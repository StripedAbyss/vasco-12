#include "HybridManufacturing.h"
#include "AccessibilityVisualization.h"
#include<stdlib.h>
#include "vasco/core/ContactTriangulation.h"
#include "vasco/core/ContourContainment.h"
#include "vasco/core/MeshValidation.h"
#include "vasco/core/PatchBoundaryStitching.h"
#include "vasco/core/SlicerMeshAdapter.h"
#include<igl/readOBJ.h>
#include<igl/writeOBJ.h>
#include <array>
#include <CGAL/IO/polygon_mesh_io.h>
#include <direct.h>
#include <cstdint>
#include <future>
#include <iomanip>
#include <limits>
#include <random>
#include <set>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

namespace
{
	// Removes structures registered by the previous input model while keeping the Polyscope window initialized.
	void ResetPolyscopeSceneForGlobalDecomposition()
	{
		polyscope::removeAllStructures();
		std::cout << "[HybridManufacturing::subtractive_accessibility_decomposition_global] "
			<< "Cleared previous Polyscope structures." << std::endl;
	}

	inline void PrepareToolForCollision(cutter& tool)
	{
		tool.cylinder_height_threshold = tool.cylinder_height + tool.ball_r;
		tool.carriage_check_radius_sq = (tool.carriage_r + 5.0) * (tool.carriage_r + 5.0);
		tool.cylinder_check_radius_sq = (tool.cylinder_r + 5.0) * (tool.cylinder_r + 5.0);
		tool.cylinder_r_sq = tool.cylinder_r * tool.cylinder_r;
		tool.carriage_r_sq = tool.carriage_r * tool.carriage_r;
		tool.total_height = tool.cylinder_height + tool.ball_r + tool.carriage_height;
	}

	// Computes a normalized triangle-face normal for the self-support diagnostics.
	bool ComputeSelfSupportFaceNormal(
		const SurfaceMesh& mesh,
		SurfaceMesh::Face_index face_index,
		Eigen::Vector3d& normal)
	{
		if (face_index == SurfaceMesh::null_face()
			|| static_cast<std::size_t>(face_index.idx()) >= mesh.number_of_faces()) {
			return false;
		}

		const auto h0 = mesh.halfedge(face_index);
		const auto h1 = mesh.next(h0);
		const auto h2 = mesh.next(h1);
		const Point_3& p0 = mesh.point(mesh.target(h0));
		const Point_3& p1 = mesh.point(mesh.target(h1));
		const Point_3& p2 = mesh.point(mesh.target(h2));
		const Eigen::Vector3d e1(
			CGAL::to_double(p1.x() - p0.x()),
			CGAL::to_double(p1.y() - p0.y()),
			CGAL::to_double(p1.z() - p0.z()));
		const Eigen::Vector3d e2(
			CGAL::to_double(p2.x() - p0.x()),
			CGAL::to_double(p2.y() - p0.y()),
			CGAL::to_double(p2.z() - p0.z()));
		normal = e1.cross(e2);
		if (!normal.allFinite() || normal.squaredNorm() <= 1e-24) {
			return false;
		}
		normal.normalize();
		return true;
	}

	// Counts per-layer downward normals from the same stored array used by the search decision.
	void PopulateLayerSelfSupportDiagnostics(
		const std::vector<SurfaceMeshSliceData>& slices,
		Layer_Graph& layer_graph)
	{
		const int node_count = layer_graph.total_node_num;
		layer_graph.self_support_associated_face_count.assign(node_count, 0);
		layer_graph.self_support_stored_normal_count.assign(node_count, 0);
		layer_graph.self_support_valid_normal_count.assign(node_count, 0);
		layer_graph.self_support_invalid_face_count.assign(node_count, 0);
		layer_graph.self_support_too_downward_face_count.assign(node_count, 0);
		layer_graph.self_support_too_downward_face_ids.assign(node_count, {});

		const Eigen::Vector3d base_normal(0.0, 0.0, 1.0);
		const double normal_z_threshold = -std::sin(PI / 3.6);
		int stored_normal_count = 0;
		const std::vector<Eigen::Vector3d>* stored_normals = nullptr;
		for (const auto& slice : slices) {
			stored_normal_count += static_cast<int>(slice.face_normals.size());
			if (stored_normals == nullptr && !slice.face_normals.empty()) {
				stored_normals = &slice.face_normals;
			}
		}
		for (int node_id = 0; node_id < node_count; ++node_id) {
			const auto index_it = layer_graph.data.index.find(node_id);
			if (index_it == layer_graph.data.index.end()) {
				continue;
			}

			const int slice_id = index_it->second.first;
			const int component_id = index_it->second.second;
			if (slice_id < 0
				|| slice_id >= static_cast<int>(slices.size())
				|| slice_id >= static_cast<int>(layer_graph.data.source_contour_ids.size())
				|| slice_id >= static_cast<int>(layer_graph.data.source_hole_contour_ids.size())
				|| component_id < 0
				|| component_id >= static_cast<int>(layer_graph.data.source_contour_ids[slice_id].size())
				|| component_id >= static_cast<int>(layer_graph.data.source_hole_contour_ids[slice_id].size())) {
				continue;
			}

			const auto& slice = slices[slice_id];
			layer_graph.self_support_stored_normal_count[node_id] =
				stored_normal_count;
			auto process_source_contour = [&](int source_contour_id) {
				if (source_contour_id < 0
					|| source_contour_id >= static_cast<int>(slice.contour_face_ids.size())) {
					return;
				}
				for (SurfaceMesh::Face_index face_index : slice.contour_face_ids[source_contour_id]) {
					++layer_graph.self_support_associated_face_count[node_id];
					const int normal_index = static_cast<int>(face_index.idx());
					if (stored_normals == nullptr
						|| normal_index < 0
						|| normal_index >= static_cast<int>(stored_normals->size())) {
						++layer_graph.self_support_invalid_face_count[node_id];
						continue;
					}
					const Eigen::Vector3d& face_normal = (*stored_normals)[normal_index];
					if (!face_normal.allFinite() || face_normal.squaredNorm() <= 1e-24) {
						++layer_graph.self_support_invalid_face_count[node_id];
						continue;
					}
					++layer_graph.self_support_valid_normal_count[node_id];
					if (face_normal.dot(base_normal) < normal_z_threshold) {
						++layer_graph.self_support_too_downward_face_count[node_id];
						layer_graph.self_support_too_downward_face_ids[node_id].push_back(
							normal_index);
					}
				}
			};

			process_source_contour(
				layer_graph.data.source_contour_ids[slice_id][component_id]);
			for (int hole_source_id :
				layer_graph.data.source_hole_contour_ids[slice_id][component_id]) {
				process_source_contour(hole_source_id);
			}
		}
	}

	inline Eigen::Vector3d ToVector3(const Eigen::MatrixXd& vec)
	{
		return { vec(0, 0), vec(1, 0), vec(2, 0) };
	}

	inline Eigen::Vector3d ComputeToolCenter(const Eigen::Vector3d& point, const Eigen::Vector3d& normal, double radius)
	{
		return point + radius * normal;
	}

	void EnsureDirectory(const std::string& directory)
	{
		if (!directory.empty()) {
			_mkdir(directory.c_str());
		}
	}

	void EnsureParentDirectory(const std::string& file_path)
	{
		const std::size_t slash_pos = file_path.find_last_of("/\\");
		if (slash_pos != std::string::npos) {
			EnsureDirectory(file_path.substr(0, slash_pos));
		}
	}

	std::string JoinPath(const std::string& directory, const std::string& file_name)
	{
		if (directory.empty()) {
			return file_name;
		}
		const char last = directory.back();
		if (last == '\\' || last == '/') {
			return directory + file_name;
		}
		return directory + "\\" + file_name;
	}

	std::string TrimCopy(const std::string& text)
	{
		const std::size_t first = text.find_first_not_of(" \t\r\n");
		if (first == std::string::npos) {
			return "";
		}
		const std::size_t last = text.find_last_not_of(" \t\r\n");
		return text.substr(first, last - first + 1);
	}

	std::string ExtractReportValue(const std::string& line, const std::string& key)
	{
		const std::size_t key_pos = line.find(key);
		if (key_pos == std::string::npos) {
			return "";
		}
		const std::size_t value_begin = key_pos + key.size();
		const std::size_t value_end = line.find(", ", value_begin);
		if (value_end == std::string::npos) {
			return TrimCopy(line.substr(value_begin));
		}
		return TrimCopy(line.substr(value_begin, value_end - value_begin));
	}

	std::string FileStemFromPath(const std::string& path)
	{
		const std::size_t slash_pos = path.find_last_of("/\\");
		std::string name = (slash_pos == std::string::npos) ? path : path.substr(slash_pos + 1);
		const std::size_t dot_pos = name.find_last_of('.');
		if (dot_pos != std::string::npos) {
			name = name.substr(0, dot_pos);
		}
		return name;
	}

	std::string TryDerivePatchFileFromSourceInput(
		const std::string& source_input_file,
		const std::string& model_file_name,
		const std::string& patch_dir)
	{
		const std::string source_stem = FileStemFromPath(source_input_file);
		const std::string model_stem = FileStemFromPath(model_file_name);
		const std::string prefix = model_stem + "-";
		if (source_stem.rfind(prefix, 0) != 0) {
			return "";
		}

		std::string token = source_stem.substr(prefix.size());
		if (token == "0_0") {
			return "";
		}

		std::vector<std::string> parts;
		std::size_t start = 0;
		while (start <= token.size()) {
			const std::size_t pos = token.find('_', start);
			if (pos == std::string::npos) {
				parts.push_back(token.substr(start));
				break;
			}
			parts.push_back(token.substr(start, pos - start));
			start = pos + 1;
		}
		if (parts.size() < 2) {
			return "";
		}

		return JoinPath(patch_dir, "block_patch-" + parts[0] + "_" + parts[1] + ".obj");
	}

	std::vector<std::string> LoadFirstFinalNodeAncestorPatchFiles(
		const std::string& ancestor_source_report_file,
		const std::string& model_file_name,
		const std::string& patch_dir)
	{
		std::ifstream report(ancestor_source_report_file);
		if (!report.is_open()) {
			std::cout << "[Warn] cannot open ancestor source report: "
				<< ancestor_source_report_file << std::endl;
			return {};
		}

		std::vector<std::string> patch_files_depth_from_final;
		std::vector<std::string> source_files_depth_from_final;
		bool in_first_final_node = false;
		std::string line;
		while (std::getline(report, line)) {
			if (line.rfind("final_node:", 0) == 0) {
				if (in_first_final_node) {
					break;
				}
				in_first_final_node = true;
				continue;
			}
			if (!in_first_final_node) {
				continue;
			}

			const std::string patch_file = ExtractReportValue(line, "patch_output_file:");
			if (!patch_file.empty() && patch_file != "<none>") {
				patch_files_depth_from_final.push_back(patch_file);
			}

			const std::string source_file = ExtractReportValue(line, "source_input_file:");
			if (!source_file.empty()) {
				source_files_depth_from_final.push_back(source_file);
			}
		}

		if (patch_files_depth_from_final.empty()) {
			for (const auto& source_file : source_files_depth_from_final) {
				const std::string patch_file =
					TryDerivePatchFileFromSourceInput(source_file, model_file_name, patch_dir);
				if (!patch_file.empty()) {
					patch_files_depth_from_final.push_back(patch_file);
				}
			}
			if (!patch_files_depth_from_final.empty()) {
				std::cout << "[Warn] ancestor report has no patch_output_file entries; "
					<< "derived patch files from source_input_file as a fallback." << std::endl;
			}
		}

		std::reverse(patch_files_depth_from_final.begin(), patch_files_depth_from_final.end());
		patch_files_depth_from_final.erase(
			std::unique(patch_files_depth_from_final.begin(), patch_files_depth_from_final.end()),
			patch_files_depth_from_final.end());
		return patch_files_depth_from_final;
	}

	// Removes vertices that are not referenced by any patch face while preserving
	// face order, vertex order, and all valid face geometry.
	bool CompactPatchToReferencedVertices(Slicer_2& patch, std::size_t& removed_vertex_count)
	{
		removed_vertex_count = 0;
		std::vector<bool> referenced_vertices(patch.positions.size(), false);
		for (const auto& triangle : patch.triangles) {
			for (int corner = 0; corner < 3; ++corner) {
				const int vertex_index = triangle[corner];
				if (vertex_index < 0 || vertex_index >= static_cast<int>(patch.positions.size())) {
					return false;
				}
				referenced_vertices[vertex_index] = true;
			}
		}

		std::vector<int> old_to_new(patch.positions.size(), -1);
		std::vector<Slicer_2::Vec3> compact_positions;
		compact_positions.reserve(patch.positions.size());
		for (std::size_t old_index = 0; old_index < patch.positions.size(); ++old_index) {
			if (!referenced_vertices[old_index]) {
				continue;
			}
			old_to_new[old_index] = static_cast<int>(compact_positions.size());
			compact_positions.push_back(patch.positions[old_index]);
		}

		for (auto& triangle : patch.triangles) {
			for (int corner = 0; corner < 3; ++corner) {
				triangle[corner] = old_to_new[triangle[corner]];
			}
		}

		removed_vertex_count = patch.positions.size() - compact_positions.size();
		patch.positions.swap(compact_positions);
		return true;
	}

	// Resolves cross-patch boundary T-junctions before graph-cut adjacency is built
	// and reports the topology changes made by the stitching module.
	bool StitchMergedPatchBoundariesForGraphCut(
		Slicer_2& merged_patch,
		std::vector<int>& merged_face_source_patch_id,
		const std::string& context)
	{
		vasco::patch_boundary_stitching::StitchOptions options;
		vasco::patch_boundary_stitching::StitchStats stats;
		const bool succeeded = vasco::patch_boundary_stitching::StitchPatchBoundariesWithTolerance(
			merged_patch,
			merged_face_source_patch_id,
			options,
			stats);

		std::cout << "[PatchBoundaryStitching] " << context
			<< ": success=" << succeeded
			<< ", model_diagonal=" << stats.model_diagonal
			<< ", model_tolerance=" << stats.model_tolerance
			<< ", initial_boundary_edges=" << stats.initial_boundary_edge_count
			<< ", initial_t_junctions=" << stats.initial_t_junction_count
			<< ", projected_vertices=" << stats.projected_vertex_count
			<< ", split_edges=" << stats.split_edge_count
			<< ", split_faces=" << stats.split_face_count
			<< ", added_faces=" << stats.added_face_count
			<< ", rejected_candidates=" << stats.rejected_candidate_count
			<< ", remaining_t_junctions=" << stats.remaining_t_junction_count
			<< ", final_boundary_edges=" << stats.final_boundary_edge_count
			<< ", iterations=" << stats.iteration_count
			<< std::endl;

		if (!succeeded || stats.remaining_t_junction_count != 0) {
			std::cerr << "[PatchBoundaryStitching] Refuse to build graph-cut adjacency from "
				<< "an invalid or incompletely stitched mesh: " << context << std::endl;
			return false;
		}
		return true;
	}

	struct MergedPatchGraphCutEdgeAnnotation {
		int first_vertex = -1;
		int second_vertex = -1;
		int first_face = -1;
		int second_face = -1;
		int incident_face_count = 0;
		int first_source_block = -1;
		int second_source_block = -1;
		double length = 0.0;
		bool is_additive_block_boundary = false;
		int base_smooth_cost = 0;
		int block_boundary_overlap_cost = 0;
	};

	struct MergedPatchGraphCutAdjacency {
		std::vector<std::vector<int>> face_relations;
		// GeneralGraph_DArraySArraySpatVarying divides these values by 10.
		std::vector<int> edge_weights_times_ten;
		std::vector<MergedPatchGraphCutEdgeAnnotation> edge_annotations;
		std::size_t additive_block_boundary_edge_count = 0;
		double additive_block_boundary_total_length = 0.0;
		std::size_t exterior_edge_count = 0;
		std::size_t non_manifold_edge_count = 0;
	};

	// Builds the shared-edge graph once for all merged-patch graph-cut modes.
	// An internal edge is an additive-block boundary when its two source block ids differ.
	bool BuildMergedPatchGraphCutAdjacency(
		const Slicer_2& merged_patch,
		const std::vector<int>& face_source_block_ids,
		double block_boundary_overlap_cost_per_unit,
		MergedPatchGraphCutAdjacency& result)
	{
		result = {};
		if (merged_patch.triangles.size() != face_source_block_ids.size()
			|| !std::isfinite(block_boundary_overlap_cost_per_unit)
			|| block_boundary_overlap_cost_per_unit < 0.0) {
			std::cerr << "[GraphCut] Cannot build merged-patch adjacency: invalid source-id count "
				<< "or block-boundary overlap weight." << std::endl;
			return false;
		}

		struct EdgeUse {
			int first_vertex = -1;
			int second_vertex = -1;
			std::array<int, 2> faces = { -1, -1 };
			int incident_face_count = 0;
		};

		auto make_edge_key = [](int first_vertex, int second_vertex) {
			const std::uint32_t low = static_cast<std::uint32_t>(
				std::min(first_vertex, second_vertex));
			const std::uint32_t high = static_cast<std::uint32_t>(
				std::max(first_vertex, second_vertex));
			return (static_cast<std::uint64_t>(low) << 32)
				| static_cast<std::uint64_t>(high);
			};

		std::unordered_map<std::uint64_t, EdgeUse> edge_uses;
		edge_uses.reserve(merged_patch.triangles.size() * 2 + 1);
		for (int face_id = 0;
			face_id < static_cast<int>(merged_patch.triangles.size());
			++face_id) {
			const auto& triangle = merged_patch.triangles[face_id];
			for (int corner = 0; corner < 3; ++corner) {
				const int first_vertex = triangle[corner];
				const int second_vertex = triangle[(corner + 1) % 3];
				if (first_vertex < 0
					|| second_vertex < 0
					|| first_vertex >= static_cast<int>(merged_patch.positions.size())
					|| second_vertex >= static_cast<int>(merged_patch.positions.size())
					|| first_vertex == second_vertex) {
					std::cerr << "[GraphCut] Cannot build merged-patch adjacency: invalid edge in face "
						<< face_id << "." << std::endl;
					return false;
				}

				const std::uint64_t key =
					make_edge_key(first_vertex, second_vertex);
				auto insertion = edge_uses.try_emplace(key);
				EdgeUse& edge = insertion.first->second;
				if (insertion.second) {
					edge.first_vertex = std::min(first_vertex, second_vertex);
					edge.second_vertex = std::max(first_vertex, second_vertex);
				}
				if (edge.incident_face_count < 2) {
					edge.faces[edge.incident_face_count] = face_id;
				}
				++edge.incident_face_count;
			}
		}

		std::vector<EdgeUse> sorted_edges;
		sorted_edges.reserve(edge_uses.size());
		for (const auto& edge_entry : edge_uses) {
			sorted_edges.push_back(edge_entry.second);
		}
		std::sort(
			sorted_edges.begin(),
			sorted_edges.end(),
			[](const EdgeUse& lhs, const EdgeUse& rhs) {
				return std::tie(lhs.first_vertex, lhs.second_vertex)
					< std::tie(rhs.first_vertex, rhs.second_vertex);
			});

		result.edge_annotations.reserve(sorted_edges.size());
		result.face_relations.reserve(sorted_edges.size());
		result.edge_weights_times_ten.reserve(sorted_edges.size());
		for (const EdgeUse& edge : sorted_edges) {
			MergedPatchGraphCutEdgeAnnotation annotation;
			annotation.first_vertex = edge.first_vertex;
			annotation.second_vertex = edge.second_vertex;
			annotation.first_face = edge.faces[0];
			annotation.second_face = edge.faces[1];
			annotation.incident_face_count = edge.incident_face_count;
			if (annotation.first_face >= 0) {
				annotation.first_source_block =
					face_source_block_ids[annotation.first_face];
			}
			if (annotation.second_face >= 0) {
				annotation.second_source_block =
					face_source_block_ids[annotation.second_face];
			}

			const auto& first_point =
				merged_patch.positions[annotation.first_vertex];
			const auto& second_point =
				merged_patch.positions[annotation.second_vertex];
			const double dx = first_point[0] - second_point[0];
			const double dy = first_point[1] - second_point[1];
			const double dz = first_point[2] - second_point[2];
			annotation.length = std::sqrt(dx * dx + dy * dy + dz * dz);
			if (!std::isfinite(annotation.length)) {
				std::cerr << "[GraphCut] Cannot build merged-patch adjacency: non-finite edge length."
					<< std::endl;
				return false;
			}

			if (edge.incident_face_count == 1) {
				++result.exterior_edge_count;
			}
			else if (edge.incident_face_count > 2) {
				++result.non_manifold_edge_count;
			}
			else {
				annotation.is_additive_block_boundary =
					annotation.first_source_block != annotation.second_source_block;
				if (annotation.is_additive_block_boundary) {
					++result.additive_block_boundary_edge_count;
					result.additive_block_boundary_total_length += annotation.length;
				}

				const double scaled_base_cost = annotation.length * 100.0;
				const int base_weight_times_ten = static_cast<int>(std::min(
					scaled_base_cost,
					static_cast<double>(std::numeric_limits<int>::max())));
				annotation.base_smooth_cost =
					std::max(1, base_weight_times_ten) / 10;

				if (annotation.is_additive_block_boundary) {
					const double overlap_cost =
						annotation.length * block_boundary_overlap_cost_per_unit;
					annotation.block_boundary_overlap_cost =
						static_cast<int>(std::min(
							static_cast<double>(std::numeric_limits<int>::max()),
							std::floor(overlap_cost + 0.5)));
				}

				const long long total_pairwise_cost =
					static_cast<long long>(annotation.base_smooth_cost)
					+ annotation.block_boundary_overlap_cost;
				const int total_weight_times_ten = static_cast<int>(std::min(
					total_pairwise_cost * 10,
					static_cast<long long>(std::numeric_limits<int>::max())));
				result.face_relations.push_back({
					annotation.first_face,
					annotation.second_face
					});
				result.edge_weights_times_ten.push_back(total_weight_times_ten);
			}
			result.edge_annotations.push_back(annotation);
		}

		std::cout << "[GraphCut] Merged-patch edge classification: edges="
			<< result.edge_annotations.size()
			<< ", adjacency_edges=" << result.face_relations.size()
			<< ", additive_block_boundary_edges="
			<< result.additive_block_boundary_edge_count
			<< ", additive_block_boundary_length="
			<< result.additive_block_boundary_total_length
			<< ", exterior_edges=" << result.exterior_edge_count
			<< ", non_manifold_edges=" << result.non_manifold_edge_count
			<< ", overlap_cost_per_unit="
			<< block_boundary_overlap_cost_per_unit
			<< std::endl;
		return result.non_manifold_edge_count == 0;
	}

	// Writes every merged-mesh edge together with its additive-block-boundary flag
	// and the overlap cost contributed by the final graph-cut labeling.
	void WriteMergedPatchGraphCutEdgeCostReport(
		const std::string& output_file,
		const std::string& summary_file,
		const MergedPatchGraphCutAdjacency& adjacency,
		const std::vector<int>& labels)
	{
		EnsureParentDirectory(output_file);
		std::ofstream report(output_file);
		const bool write_edge_report = report.is_open();
		if (!report.is_open()) {
			std::cerr << "[GraphCut] Cannot write edge cost report: "
				<< output_file << std::endl;
		}

		if (write_edge_report) {
			report << std::setprecision(17);
			report
				<< "edge_id,vertex_0,vertex_1,incident_face_count,"
				<< "face_0,face_1,source_block_0,source_block_1,"
				<< "is_additive_block_boundary,length,"
				<< "label_0,label_1,is_graphcut_patch_boundary,"
				<< "base_smooth_cost,block_boundary_overlap_cost,"
				<< "overlap_cost_contribution\n";
		}

		double overlapping_boundary_length = 0.0;
		long long overlap_cost = 0;
		std::size_t overlapping_edge_count = 0;
		for (std::size_t edge_id = 0;
			edge_id < adjacency.edge_annotations.size();
			++edge_id) {
			const auto& edge = adjacency.edge_annotations[edge_id];
			const int first_label =
				edge.first_face >= 0
				&& edge.first_face < static_cast<int>(labels.size())
				? labels[edge.first_face]
				: -1;
			const int second_label =
				edge.second_face >= 0
				&& edge.second_face < static_cast<int>(labels.size())
				? labels[edge.second_face]
				: -1;
			const bool is_graphcut_patch_boundary =
				edge.incident_face_count == 2
				&& first_label >= 0
				&& second_label >= 0
				&& first_label != second_label;
			const int contribution =
				edge.is_additive_block_boundary && is_graphcut_patch_boundary
				? edge.block_boundary_overlap_cost
				: 0;
			if (edge.is_additive_block_boundary && is_graphcut_patch_boundary) {
				++overlapping_edge_count;
				overlapping_boundary_length += edge.length;
				overlap_cost += contribution;
			}

			if (write_edge_report) {
				report
					<< edge_id << ","
					<< edge.first_vertex << ","
					<< edge.second_vertex << ","
					<< edge.incident_face_count << ","
					<< edge.first_face << ","
					<< edge.second_face << ","
					<< edge.first_source_block << ","
					<< edge.second_source_block << ","
					<< (edge.is_additive_block_boundary ? 1 : 0) << ","
					<< edge.length << ","
					<< first_label << ","
					<< second_label << ","
					<< (is_graphcut_patch_boundary ? 1 : 0) << ","
					<< edge.base_smooth_cost << ","
					<< edge.block_boundary_overlap_cost << ","
					<< contribution << "\n";
			}
		}

		std::cout << "[GraphCut] Additive-block/graph-cut boundary overlap: edges="
			<< overlapping_edge_count
			<< ", length=" << overlapping_boundary_length
			<< ", cost=" << overlap_cost
			<< ", report=" << output_file
			<< std::endl;

		EnsureParentDirectory(summary_file);
		std::ofstream summary(summary_file);
		if (!summary.is_open()) {
			std::cerr << "[GraphCut] Cannot write boundary overlap summary: "
				<< summary_file << std::endl;
			return;
		}
		const double overlap_ratio =
			adjacency.additive_block_boundary_total_length > 0.0
			? overlapping_boundary_length
				/ adjacency.additive_block_boundary_total_length
			: 0.0;
		summary << std::setprecision(17);
		summary << "L_block="
			<< adjacency.additive_block_boundary_total_length << "\n";
		summary << "L_overlap=" << overlapping_boundary_length << "\n";
		summary << "R_overlap=" << overlap_ratio << "\n";
		std::cout << "[GraphCut] Boundary overlap summary saved: "
			<< summary_file << std::endl;
	}

	// Makes every orientation of one source block feasible for a face while
	// keeping all labels belonging to other blocks unchanged.
	bool EnableAllOrientationLabelsForSourceBlock(
		std::vector<int>& face_data_costs,
		int source_block_id,
		int orientation_count)
	{
		if (source_block_id <= 0 || orientation_count <= 0) {
			return false;
		}

		const std::size_t first_label =
			static_cast<std::size_t>(source_block_id - 1) * static_cast<std::size_t>(orientation_count);
		const std::size_t end_label = first_label + static_cast<std::size_t>(orientation_count);
		if (end_label > face_data_costs.size()) {
			return false;
		}

		std::fill(
			face_data_costs.begin() + first_label,
			face_data_costs.begin() + end_label,
			0);
		return true;
	}

	struct InitializedGraphCutResult {
		std::vector<int> labels;
		long long initial_energy = std::numeric_limits<long long>::max();
		long long energy = std::numeric_limits<long long>::max();
		long long data_energy = 0;
		long long smooth_energy = 0;
		long long label_energy = 0;
		bool succeeded = false;
	};

	// Runs the same GCO objective as GeneralGraph_DArraySArraySpatVarying,
	// but starts alpha-expansion from a caller-provided valid labeling.
	InitializedGraphCutResult RunGeneralGraphCutWithInitialLabels(
		int node_count,
		int label_count,
		const std::vector<std::vector<int>>& data_costs,
		const std::vector<std::vector<int>>& adjacency,
		const std::vector<int>& edge_lengths,
		const std::vector<int>& initial_labels,
		const std::string& initialization_name,
		bool log_progress = true)
	{
		InitializedGraphCutResult result;
		result.labels.assign(static_cast<std::size_t>(std::max(0, node_count)), 0);
		if (node_count <= 0
			|| label_count <= 0
			|| data_costs.size() != static_cast<std::size_t>(node_count)
			|| adjacency.size() != edge_lengths.size()
			|| initial_labels.size() != static_cast<std::size_t>(node_count)) {
			std::cerr << "[GraphCut] Cannot apply the supplied initialization because the graph dimensions are invalid."
				<< std::endl;
			return result;
		}

		for (int label : initial_labels) {
			if (label < 0 || label >= label_count) {
				std::cerr << "[GraphCut] Cannot apply the supplied initialization because it contains label "
					<< label << " outside [0, " << (label_count - 1) << "]." << std::endl;
				return result;
			}
		}

		std::vector<int> packed_data(
			static_cast<std::size_t>(node_count) * static_cast<std::size_t>(label_count));
		for (int node = 0; node < node_count; ++node) {
			if (data_costs[static_cast<std::size_t>(node)].size()
				!= static_cast<std::size_t>(label_count)) {
				std::cerr << "[GraphCut] Cannot apply the supplied initialization because data-cost row "
					<< node << " has an invalid size." << std::endl;
				return result;
			}
			for (int label = 0; label < label_count; ++label) {
				packed_data[
					static_cast<std::size_t>(node) * static_cast<std::size_t>(label_count)
					+ static_cast<std::size_t>(label)] =
					data_costs[static_cast<std::size_t>(node)][static_cast<std::size_t>(label)];
			}
		}

		std::vector<int> smooth_costs(
			static_cast<std::size_t>(label_count) * static_cast<std::size_t>(label_count));
		for (int first_label = 0; first_label < label_count; ++first_label) {
			for (int second_label = 0; second_label < label_count; ++second_label) {
				smooth_costs[
					static_cast<std::size_t>(first_label)
					+ static_cast<std::size_t>(second_label) * static_cast<std::size_t>(label_count)] =
					(first_label == second_label) ? 0 : 1;
			}
		}
		std::vector<int> label_costs(static_cast<std::size_t>(label_count), 100000);

		try {
			GCoptimizationGeneralGraph graph_cut(node_count, label_count);
			graph_cut.setDataCost(packed_data.data());
			graph_cut.setSmoothCost(smooth_costs.data());
			graph_cut.setLabelCost(label_costs.data());

			for (std::size_t edge_index = 0; edge_index < adjacency.size(); ++edge_index) {
				if (adjacency[edge_index].size() != 2) {
					std::cerr << "[GraphCut] Skip malformed adjacency entry " << edge_index << "."
						<< std::endl;
					continue;
				}
				graph_cut.setNeighbors(
					adjacency[edge_index][0],
					adjacency[edge_index][1],
					edge_lengths[edge_index] / 10);
			}
			for (int node = 0; node < node_count; ++node) {
				graph_cut.setLabel(node, initial_labels[static_cast<std::size_t>(node)]);
			}

			if (log_progress) {
				std::cout << "[GraphCut] Applied " << initialization_name
					<< " as the initial labeling." << std::endl;
			}
			result.initial_energy = graph_cut.compute_energy();
			if (log_progress) {
				std::cout << "Before optimization energy is " << result.initial_energy << std::endl;
			}
			graph_cut.expansion(20);
			result.energy = graph_cut.compute_energy();
			result.data_energy = graph_cut.giveDataEnergy();
			result.smooth_energy = graph_cut.giveSmoothEnergy();
			result.label_energy = graph_cut.giveLabelEnergy();
			if (log_progress) {
				std::cout << "After optimization energy is " << result.energy << std::endl;
				std::cout << "***********" << std::endl;
				std::cout << "DataEnergy: " << result.data_energy << std::endl;
				std::cout << "SmoothEnergy: " << result.smooth_energy << std::endl;
				std::cout << "LabelEnergy: " << result.label_energy << std::endl;
				std::cout << "***********" << std::endl;
			}

			for (int node = 0; node < node_count; ++node) {
				result.labels[static_cast<std::size_t>(node)] = graph_cut.whatLabel(node);
			}
			result.succeeded = true;
		}
		catch (GCException& exception) {
			exception.Report();
		}

		return result;
	}

	struct RandomRestartSummary {
		int restart = 0;
		unsigned int seed = 0;
		bool succeeded = false;
		long long initial_energy = 0;
		long long energy = 0;
		long long data_energy = 0;
		long long smooth_energy = 0;
		long long label_energy = 0;
		int label_count = 0;
	};

	struct RandomRestartTaskResult {
		int restart = 0;
		unsigned int seed = 0;
		InitializedGraphCutResult graph_cut_result;
	};

	// Derives a stable, distinct seed for one restart from the user-provided
	// base seed and restart index.
	unsigned int DeriveRandomRestartSeed(unsigned int base_seed, int restart)
	{
		std::uint32_t value =
			static_cast<std::uint32_t>(base_seed)
			+ 0x9e3779b9u * (static_cast<std::uint32_t>(restart) + 1u);
		value ^= value >> 16;
		value *= 0x85ebca6bu;
		value ^= value >> 13;
		value *= 0xc2b2ae35u;
		value ^= value >> 16;
		return static_cast<unsigned int>(value);
	}

	// Resolves the requested worker count. Automatic mode is capped at four
	// because each concurrent GCO instance owns sizeable graph work buffers.
	int ResolveRandomRestartThreadCount(int requested_thread_count, int restart_count)
	{
		if (restart_count <= 0) {
			return 1;
		}
		if (requested_thread_count > 0) {
			return std::max(1, std::min(requested_thread_count, restart_count));
		}

		const unsigned int hardware_threads = std::thread::hardware_concurrency();
		const int automatic_thread_count =
			hardware_threads == 0 ? 1 : static_cast<int>(hardware_threads);
		return std::max(1, std::min({ automatic_thread_count, restart_count, 4 }));
	}

	// Runs graph cut repeatedly from uniformly sampled feasible labels and
	// returns the successful result with the lowest final total energy.
	InitializedGraphCutResult RunRandomFeasibleLabelGraphCutRestarts(
		int node_count,
		int label_count,
		const std::vector<std::vector<int>>& data_costs,
		const std::vector<std::vector<int>>& adjacency,
		const std::vector<int>& edge_lengths,
		int restart_count,
		unsigned int random_seed,
		int requested_thread_count,
		const std::string& context,
		const std::string& report_file)
	{
		InitializedGraphCutResult best_result;
		if (node_count <= 0
			|| label_count <= 0
			|| restart_count <= 0
			|| data_costs.size() != static_cast<std::size_t>(node_count)) {
			std::cerr << "[GraphCut] Cannot run random feasible-label restarts for "
				<< context << " because the graph dimensions are invalid." << std::endl;
			return best_result;
		}

		// Build each node's feasible-label list once. This avoids rescanning all
		// labels during every random restart.
		std::vector<std::vector<int>> feasible_labels(static_cast<std::size_t>(node_count));
		for (int node = 0; node < node_count; ++node) {
			const auto& node_costs = data_costs[static_cast<std::size_t>(node)];
			if (node_costs.size() != static_cast<std::size_t>(label_count)) {
				std::cerr << "[GraphCut] Cannot run random feasible-label restarts for "
					<< context << " because data-cost row " << node
					<< " has an invalid size." << std::endl;
				return best_result;
			}
			auto& node_feasible_labels = feasible_labels[static_cast<std::size_t>(node)];
			for (int label = 0; label < label_count; ++label) {
				if (node_costs[static_cast<std::size_t>(label)] == 0) {
					node_feasible_labels.push_back(label);
				}
			}
			if (node_feasible_labels.empty()) {
				std::cerr << "[GraphCut] Cannot run random feasible-label restarts for "
					<< context << " because node " << node
					<< " has no feasible label." << std::endl;
				return best_result;
			}
		}

		std::vector<RandomRestartSummary> restart_summaries;
		restart_summaries.reserve(static_cast<std::size_t>(restart_count));
		int best_restart = -1;
		const int thread_count =
			ResolveRandomRestartThreadCount(requested_thread_count, restart_count);
		gcoclock(); // Initialize GCO's Windows timer frequency before workers start.
		std::cout << "[GraphCut] Run " << restart_count << ' ' << context
			<< " random restarts with " << thread_count
			<< " worker thread(s), base_seed=" << random_seed << '.' << std::endl;

		for (int batch_begin = 0; batch_begin < restart_count; batch_begin += thread_count) {
			const int batch_end = std::min(restart_count, batch_begin + thread_count);
			std::vector<std::future<RandomRestartTaskResult>> futures;
			futures.reserve(static_cast<std::size_t>(batch_end - batch_begin));

			for (int restart = batch_begin; restart < batch_end; ++restart) {
				const unsigned int restart_seed =
					DeriveRandomRestartSeed(random_seed, restart);
				futures.push_back(std::async(
					std::launch::async,
					[&, restart, restart_seed]() {
						std::mt19937 random_engine(restart_seed);
						std::vector<int> initial_labels(
							static_cast<std::size_t>(node_count),
							0);
						for (int node = 0; node < node_count; ++node) {
							const auto& node_feasible_labels =
								feasible_labels[static_cast<std::size_t>(node)];
							std::uniform_int_distribution<std::size_t> distribution(
								0,
								node_feasible_labels.size() - 1);
							initial_labels[static_cast<std::size_t>(node)] =
								node_feasible_labels[distribution(random_engine)];
						}

						const std::string initialization_name =
							context + " random restart " + std::to_string(restart + 1)
							+ "/" + std::to_string(restart_count)
							+ " (seed=" + std::to_string(restart_seed) + ")";
						RandomRestartTaskResult task_result;
						task_result.restart = restart;
						task_result.seed = restart_seed;
						task_result.graph_cut_result =
							RunGeneralGraphCutWithInitialLabels(
								node_count,
								label_count,
								data_costs,
								adjacency,
								edge_lengths,
								initial_labels,
								initialization_name,
								false);
						return task_result;
					}));
			}

			for (auto& future : futures) {
				const RandomRestartTaskResult task_result = future.get();
				const InitializedGraphCutResult& restart_result =
					task_result.graph_cut_result;
				const int result_label_count = restart_result.succeeded
					? static_cast<int>(std::unordered_set<int>(
						restart_result.labels.begin(),
						restart_result.labels.end()).size())
					: 0;
				restart_summaries.push_back({
					task_result.restart + 1,
					task_result.seed,
					restart_result.succeeded,
					restart_result.initial_energy,
					restart_result.energy,
					restart_result.data_energy,
					restart_result.smooth_energy,
					restart_result.label_energy,
					result_label_count
					});

				std::cout << "[GraphCut] " << context
					<< " restart " << (task_result.restart + 1)
					<< "/" << restart_count
					<< ", seed=" << task_result.seed
					<< ", succeeded=" << restart_result.succeeded;
				if (restart_result.succeeded) {
					std::cout << ", initial_energy=" << restart_result.initial_energy
						<< ", final_energy=" << restart_result.energy;
				}
				std::cout << std::endl;

				if (restart_result.succeeded
					&& restart_result.energy < best_result.energy) {
					best_result = restart_result;
					best_restart = task_result.restart;
				}
			}
		}

		if (best_restart >= 0) {
			std::cout << "[GraphCut] Best " << context
				<< " random restart=" << (best_restart + 1)
				<< "/" << restart_count
				<< ", energy=" << best_result.energy
				<< ", DataEnergy=" << best_result.data_energy
				<< ", SmoothEnergy=" << best_result.smooth_energy
				<< ", LabelEnergy=" << best_result.label_energy
				<< std::endl;
		}

		std::ofstream report(report_file);
		if (report.is_open()) {
			report
				<< "restart,restart_seed,succeeded,initial_energy,final_energy,"
				<< "data_energy,smooth_energy,label_energy,label_count,is_best\n";
			for (const auto& summary : restart_summaries) {
				report
					<< summary.restart << ','
					<< summary.seed << ','
					<< (summary.succeeded ? 1 : 0) << ','
					<< summary.initial_energy << ','
					<< summary.energy << ','
					<< summary.data_energy << ','
					<< summary.smooth_energy << ','
					<< summary.label_energy << ','
					<< summary.label_count << ','
					<< ((summary.restart == best_restart + 1) ? 1 : 0)
					<< '\n';
			}
			std::cout << "[GraphCut] Random-restart report saved to "
				<< report_file << std::endl;
		}
		else {
			std::cerr << "[GraphCut] Cannot write random-restart report: "
				<< report_file << std::endl;
		}

		return best_result;
	}

	struct LayerContainmentPolygon {
		std::vector<Point_2> points;
		Polygon_2 polygon;
		bool is_simple = false;
		bool repaired = false;
	};

	double SquaredDistance2D(const Point_2& a, const Point_2& b)
	{
		const double dx = CGAL::to_double(a.x() - b.x());
		const double dy = CGAL::to_double(a.y() - b.y());
		return dx * dx + dy * dy;
	}

	double SquaredDistancePointToSegment2D(const Point_2& p, const Point_2& a, const Point_2& b)
	{
		const double px = CGAL::to_double(p.x());
		const double py = CGAL::to_double(p.y());
		const double ax = CGAL::to_double(a.x());
		const double ay = CGAL::to_double(a.y());
		const double bx = CGAL::to_double(b.x());
		const double by = CGAL::to_double(b.y());
		const double abx = bx - ax;
		const double aby = by - ay;
		const double len_sq = abx * abx + aby * aby;
		if (len_sq <= 1e-24) {
			const double dx = px - ax;
			const double dy = py - ay;
			return dx * dx + dy * dy;
		}

		const double t = std::max(0.0, std::min(1.0, ((px - ax) * abx + (py - ay) * aby) / len_sq));
		const double cx = ax + t * abx;
		const double cy = ay + t * aby;
		const double dx = px - cx;
		const double dy = py - cy;
		return dx * dx + dy * dy;
	}

	bool IsFinitePoint2(const Point_2& point)
	{
		return std::isfinite(CGAL::to_double(point.x()))
			&& std::isfinite(CGAL::to_double(point.y()));
	}

	bool IsPointNearPolyline2D(const std::vector<Point_2>& points, const Point_2& point, double distance)
	{
		if (points.size() < 2) {
			return false;
		}

		const double distance_sq = distance * distance;
		for (std::size_t i = 0; i < points.size(); ++i) {
			const Point_2& a = points[i];
			const Point_2& b = points[(i + 1) % points.size()];
			if (SquaredDistancePointToSegment2D(point, a, b) <= distance_sq) {
				return true;
			}
		}
		return false;
	}

	bool IsPointInsideEvenOddPolygon2D(const std::vector<Point_2>& points, const Point_2& point)
	{
		if (points.size() < 3) {
			return false;
		}

		const double x = CGAL::to_double(point.x());
		const double y = CGAL::to_double(point.y());
		bool inside = false;
		for (std::size_t i = 0, j = points.size() - 1; i < points.size(); j = i++) {
			const double xi = CGAL::to_double(points[i].x());
			const double yi = CGAL::to_double(points[i].y());
			const double xj = CGAL::to_double(points[j].x());
			const double yj = CGAL::to_double(points[j].y());
			const bool crosses = ((yi > y) != (yj > y))
				&& (x < (xj - xi) * (y - yi) / (yj - yi + 1e-30) + xi);
			if (crosses) {
				inside = !inside;
			}
		}
		return inside;
	}

	bool IsNearlyCollinear2D(const Point_2& a, const Point_2& b, const Point_2& c, double eps)
	{
		return SquaredDistancePointToSegment2D(b, a, c) <= eps * eps;
	}

	std::vector<Point_2> CleanLayerPolygonPoints(const std::vector<Point_2>& raw_points, double eps)
	{
		std::vector<Point_2> points;
		points.reserve(raw_points.size());
		const double eps_sq = eps * eps;
		for (const auto& point : raw_points) {
			if (!IsFinitePoint2(point)) {
				continue;
			}
			if (!points.empty() && SquaredDistance2D(points.back(), point) <= eps_sq) {
				continue;
			}
			points.push_back(point);
		}

		if (points.size() >= 2 && SquaredDistance2D(points.front(), points.back()) <= eps_sq) {
			points.pop_back();
		}

		bool changed = true;
		while (changed && points.size() >= 3) {
			changed = false;
			std::vector<Point_2> filtered;
			filtered.reserve(points.size());
			for (std::size_t i = 0; i < points.size(); ++i) {
				const Point_2& prev = points[(i + points.size() - 1) % points.size()];
				const Point_2& current = points[i];
				const Point_2& next = points[(i + 1) % points.size()];
				if (IsNearlyCollinear2D(prev, current, next, eps)) {
					changed = true;
					continue;
				}
				filtered.push_back(current);
			}
			points.swap(filtered);
		}

		return points;
	}

	LayerContainmentPolygon PrepareLayerContainmentPolygon(
		const std::vector<Point_2>& raw_points,
		double repair_eps)
	{
		LayerContainmentPolygon result;
		result.points = CleanLayerPolygonPoints(raw_points, repair_eps);
		result.repaired = result.points.size() != raw_points.size();
		if (result.points.size() >= 3) {
			result.polygon = Polygon_2(result.points.begin(), result.points.end());
			result.is_simple = result.polygon.is_simple();
		}
		return result;
	}

	bool IsPointInsidePreparedLayerPolygon(
		const LayerContainmentPolygon& layer_polygon,
		const Point_2& point,
		double offset)
	{
		if (layer_polygon.points.size() < 3) {
			return false;
		}

		if (layer_polygon.is_simple
			&& layer_polygon.polygon.bounded_side(point) != CGAL::ON_UNBOUNDED_SIDE) {
			return true;
		}

		if (!layer_polygon.is_simple && IsPointInsideEvenOddPolygon2D(layer_polygon.points, point)) {
			return true;
		}

		return IsPointNearPolyline2D(layer_polygon.points, point, offset);
	}

	// Tests one layer component with even-odd semantics: outer minus all direct holes.
	bool IsPointInsidePreparedMaterialRegion(
		const LayerContainmentPolygon& outer_polygon,
		const std::vector<LayerContainmentPolygon>& hole_polygons,
		const Point_2& point,
		double offset)
	{
		if (!IsPointInsidePreparedLayerPolygon(outer_polygon, point, offset)) {
			return false;
		}
		for (const LayerContainmentPolygon& hole : hole_polygons) {
			if (IsPointInsidePreparedLayerPolygon(hole, point, offset)) {
				return false;
			}
		}
		return true;
	}

	bool IsPointWithinLayerZBand(double layer_z, double point_z, double layer_height, double slack)
	{
		const double z_diff = layer_z - point_z;
		return z_diff > -(layer_height + slack)
			&& z_diff <= (layer_height + slack);
	}

	// Tests the material slab generated by extruding a layer upward by one layer height.
	// A small negative tolerance keeps points numerically on the lower plane in the slab.
	bool IsPointWithinExtrudedLayerZSlab(double layer_z, double point_z, double layer_height, double slack)
	{
		const double z_diff = point_z - layer_z;
		return z_diff >= -slack
			&& z_diff <= layer_height + slack;
	}

	// Finds the nearest lower layer whose upward extrusion contains the point.
	// Layers slightly above the point are considered only as a numerical fallback.
	int FindContainingExtrudedMaterialLayer(
		const Layer_Graph& layer_graph,
		const std::vector<LayerContainmentPolygon>& outer_polygons,
		const std::vector<std::vector<LayerContainmentPolygon>>& hole_polygons,
		const Eigen::Vector3d& point,
		double layer_height,
		double z_slack,
		double boundary_offset)
	{
		int best_layer_id = -1;
		int best_z_category = 2;
		double best_z_distance = MAX_D;
		const Point_2 point_xy(point.x(), point.y());

		for (int layer_id = 0; layer_id < layer_graph.total_node_num; ++layer_id) {
			if (layer_id >= static_cast<int>(outer_polygons.size())
				|| layer_id >= static_cast<int>(hole_polygons.size())) {
				continue;
			}

			const auto index_it = layer_graph.data.index.find(layer_id);
			if (index_it == layer_graph.data.index.end()) {
				continue;
			}
			const int slice_id = index_it->second.first;
			const int contour_id = index_it->second.second;
			if (slice_id < 0
				|| slice_id >= static_cast<int>(layer_graph.data.z_value.size())
				|| contour_id < 0
				|| contour_id >= static_cast<int>(layer_graph.data.z_value[slice_id].size())
				|| layer_graph.data.z_value[slice_id][contour_id].empty()) {
				continue;
			}

			const double layer_z = layer_graph.data.z_value[slice_id][contour_id][0];
			if (!IsPointWithinExtrudedLayerZSlab(layer_z, point.z(), layer_height, z_slack)
				|| !IsPointInsidePreparedMaterialRegion(
					outer_polygons[layer_id],
					hole_polygons[layer_id],
					point_xy,
					boundary_offset)) {
				continue;
			}

			const double z_diff = point.z() - layer_z;
			const int z_category = (z_diff >= 0.0) ? 0 : 1;
			const double z_distance = std::abs(z_diff);
			if (z_category < best_z_category
				|| (z_category == best_z_category && z_distance < best_z_distance)) {
				best_layer_id = layer_id;
				best_z_category = z_category;
				best_z_distance = z_distance;
			}
		}

		return best_layer_id;
	}

	// Finds every material layer whose upward extrusion contains the point.
	// Area-S constraints use all matches so overlapping layer components cannot
	// hide an inaccessible point by assigning it only to a different component.
	std::vector<int> FindContainingExtrudedMaterialLayers(
		const Layer_Graph& layer_graph,
		const std::vector<LayerContainmentPolygon>& outer_polygons,
		const std::vector<std::vector<LayerContainmentPolygon>>& hole_polygons,
		const Eigen::Vector3d& point,
		double layer_height,
		double z_slack,
		double boundary_offset)
	{
		std::vector<int> layer_ids;
		const Point_2 point_xy(point.x(), point.y());

		for (int layer_id = 0; layer_id < layer_graph.total_node_num; ++layer_id) {
			if (layer_id >= static_cast<int>(outer_polygons.size())
				|| layer_id >= static_cast<int>(hole_polygons.size())) {
				continue;
			}

			const auto index_it = layer_graph.data.index.find(layer_id);
			if (index_it == layer_graph.data.index.end()) {
				continue;
			}
			const int slice_id = index_it->second.first;
			const int contour_id = index_it->second.second;
			if (slice_id < 0
				|| slice_id >= static_cast<int>(layer_graph.data.z_value.size())
				|| contour_id < 0
				|| contour_id >= static_cast<int>(layer_graph.data.z_value[slice_id].size())
				|| layer_graph.data.z_value[slice_id][contour_id].empty()) {
				continue;
			}

			const double layer_z = layer_graph.data.z_value[slice_id][contour_id][0];
			if (IsPointWithinExtrudedLayerZSlab(
					layer_z,
					point.z(),
					layer_height,
					z_slack)
				&& IsPointInsidePreparedMaterialRegion(
					outer_polygons[layer_id],
					hole_polygons[layer_id],
					point_xy,
					boundary_offset)) {
				layer_ids.push_back(layer_id);
			}
		}

		return layer_ids;
	}

	void WriteLineMarkerObj(
		std::ofstream& obj,
		int& next_vertex_index,
		const Eigen::Vector3d& point,
		const std::array<double, 3>& color,
		double radius)
	{
		const int x0 = next_vertex_index;
		obj << "v " << point.x() - radius << " " << point.y() << " " << point.z() << " "
			<< color[0] << " " << color[1] << " " << color[2] << "\n";
		const int x1 = next_vertex_index + 1;
		obj << "v " << point.x() + radius << " " << point.y() << " " << point.z() << " "
			<< color[0] << " " << color[1] << " " << color[2] << "\n";
		const int y0 = next_vertex_index + 2;
		obj << "v " << point.x() << " " << point.y() - radius << " " << point.z() << " "
			<< color[0] << " " << color[1] << " " << color[2] << "\n";
		const int y1 = next_vertex_index + 3;
		obj << "v " << point.x() << " " << point.y() + radius << " " << point.z() << " "
			<< color[0] << " " << color[1] << " " << color[2] << "\n";
		const int z0 = next_vertex_index + 4;
		obj << "v " << point.x() << " " << point.y() << " " << point.z() - radius << " "
			<< color[0] << " " << color[1] << " " << color[2] << "\n";
		const int z1 = next_vertex_index + 5;
		obj << "v " << point.x() << " " << point.y() << " " << point.z() + radius << " "
			<< color[0] << " " << color[1] << " " << color[2] << "\n";
		obj << "l " << x0 << " " << x1 << "\n";
		obj << "l " << y0 << " " << y1 << "\n";
		obj << "l " << z0 << " " << z1 << "\n";
		next_vertex_index += 6;
	}

	void WriteFlatSegmentQuadObj(
		std::ofstream& obj,
		int& next_vertex_index,
		const Eigen::Vector3d& a,
		const Eigen::Vector3d& b,
		const std::array<double, 3>& color,
		double width)
	{
		const double dx = b.x() - a.x();
		const double dy = b.y() - a.y();
		const double length = std::sqrt(dx * dx + dy * dy);
		if (length <= 1e-12) {
			return;
		}

		const double half_width = width * 0.5;
		const double nx = -dy / length * half_width;
		const double ny = dx / length * half_width;

		const int v0 = next_vertex_index;
		obj << "v " << a.x() + nx << " " << a.y() + ny << " " << a.z() << " "
			<< color[0] << " " << color[1] << " " << color[2] << "\n";
		const int v1 = next_vertex_index + 1;
		obj << "v " << a.x() - nx << " " << a.y() - ny << " " << a.z() << " "
			<< color[0] << " " << color[1] << " " << color[2] << "\n";
		const int v2 = next_vertex_index + 2;
		obj << "v " << b.x() - nx << " " << b.y() - ny << " " << b.z() << " "
			<< color[0] << " " << color[1] << " " << color[2] << "\n";
		const int v3 = next_vertex_index + 3;
		obj << "v " << b.x() + nx << " " << b.y() + ny << " " << b.z() << " "
			<< color[0] << " " << color[1] << " " << color[2] << "\n";
		obj << "f " << v0 << " " << v1 << " " << v2 << "\n";
		obj << "f " << v0 << " " << v2 << " " << v3 << "\n";
		next_vertex_index += 4;
	}

	void WriteDfsLayerContainmentDebugObj(
		const std::string& output_file,
		const std::vector<LayerContainmentPolygon>& layer_polygons,
		const std::vector<double>& layer_z_values,
		const std::vector<Eigen::Vector3d>& assigned_area_s_points,
		const std::vector<Eigen::Vector3d>& unassigned_area_s_points,
		double marker_radius,
		double layer_polygon_offset,
		double layer_z_slack)
	{
		EnsureParentDirectory(output_file);
		std::ofstream obj(output_file);
		if (!obj.is_open()) {
			std::cout << "[DFSLayerContainmentDebug] cannot open file for writing: "
				<< output_file << std::endl;
			return;
		}

		obj << std::setprecision(17);
		obj << "# Green/cyan/orange/red wire loops: layer polygons used by IsPointInsidePreparedLayerPolygon\n";
		obj << "# Cyan markers: area_S assigned to some layer\n";
		obj << "# Magenta markers: area_S still not assigned to any layer\n";
		obj << "# layer_polygon_offset " << layer_polygon_offset << "\n";
		obj << "# layer_z_slack " << layer_z_slack << "\n";

		int next_vertex_index = 1;
		for (std::size_t layer_id = 0; layer_id < layer_polygons.size(); ++layer_id) {
			const auto& layer_polygon = layer_polygons[layer_id];
			if (layer_polygon.points.size() < 2 || layer_id >= layer_z_values.size()) {
				continue;
			}

			const std::array<double, 3> color = layer_polygon.is_simple
				? (layer_polygon.repaired ? std::array<double, 3>{0.0, 0.85, 1.0} : std::array<double, 3>{ 0.0, 1.0, 0.0 })
				: std::array<double, 3>{ 1.0, 0.2, 0.0 };

			obj << "g layer_polygon_" << layer_id << "\n";
			obj << "# layer_id " << layer_id
				<< " z " << layer_z_values[layer_id]
				<< " point_count " << layer_polygon.points.size()
				<< " is_simple " << layer_polygon.is_simple
				<< " repaired " << layer_polygon.repaired << "\n";

			const int first_vertex = next_vertex_index;
			for (const auto& point : layer_polygon.points) {
				obj << "v "
					<< CGAL::to_double(point.x()) << " "
					<< CGAL::to_double(point.y()) << " "
					<< layer_z_values[layer_id] << " "
					<< color[0] << " " << color[1] << " " << color[2] << "\n";
				++next_vertex_index;
			}
			obj << "l";
			for (int vertex_id = first_vertex; vertex_id < next_vertex_index; ++vertex_id) {
				obj << " " << vertex_id;
			}
			obj << " " << first_vertex << "\n";

			obj << "g layer_polygon_" << layer_id << "_visible_edges\n";
			const double edge_width = std::max(0.02, marker_radius * 0.35);
			for (std::size_t edge_id = 0; edge_id < layer_polygon.points.size(); ++edge_id) {
				const Point_2& a2 = layer_polygon.points[edge_id];
				const Point_2& b2 = layer_polygon.points[(edge_id + 1) % layer_polygon.points.size()];
				const Eigen::Vector3d a(
					CGAL::to_double(a2.x()),
					CGAL::to_double(a2.y()),
					layer_z_values[layer_id]);
				const Eigen::Vector3d b(
					CGAL::to_double(b2.x()),
					CGAL::to_double(b2.y()),
					layer_z_values[layer_id]);
				WriteFlatSegmentQuadObj(obj, next_vertex_index, a, b, color, edge_width);
			}
		}

		obj << "g assigned_area_s_points\n";
		for (const auto& point : assigned_area_s_points) {
			WriteLineMarkerObj(obj, next_vertex_index, point, { 0.05, 1.0, 1.0 }, marker_radius);
		}

		obj << "g unassigned_area_s_points\n";
		for (const auto& point : unassigned_area_s_points) {
			WriteLineMarkerObj(obj, next_vertex_index, point, { 1.0, 0.0, 0.85 }, marker_radius);
		}

		std::cout << "[DFSLayerContainmentDebug] wrote layer containment debug OBJ: "
			<< output_file
			<< ", polygons=" << layer_polygons.size()
			<< ", assigned_area_s=" << assigned_area_s_points.size()
			<< ", unassigned_area_s=" << unassigned_area_s_points.size()
			<< std::endl;
	}

	std::string MakeAccessibilityDebugNodeTag(int height_of_beam_search, int cont_number_of_queue, int now_last_node)
	{
		return "_h" + std::to_string(height_of_beam_search)
			+ "_q" + std::to_string(cont_number_of_queue)
			+ "_node" + std::to_string(now_last_node);
	}

	std::string MakeAccessibilityDebugFileToken(const std::string& path)
	{
		const std::size_t slash_pos = path.find_last_of("/\\");
		std::string token = (slash_pos == std::string::npos) ? path : path.substr(slash_pos + 1);
		const std::size_t dot_pos = token.find_last_of('.');
		if (dot_pos != std::string::npos) {
			token = token.substr(0, dot_pos);
		}

		for (char& ch : token) {
			if (ch == ' ' || ch == ':' || ch == '/' || ch == '\\' || ch == '.') {
				ch = '_';
			}
		}
		return token.empty() ? "model" : token;
	}

}

void HybridManufacturing::BuildSurfaceMeshSlices(std::vector<SurfaceMeshSliceData>& slices) const
{
	slices.clear();
	const double layer_height_value = Katana::Instance().config.get("layer_height");
	if (layer_height_value <= 0) {
		return;
	}

	// Keep the normal array indexed by SurfaceMesh::Face_index because layer
	// contours store global face IDs rather than a compact local face order.
	std::vector<Eigen::Vector3d> face_normals(
		current_node_mesh_rotated.number_of_faces(),
		Eigen::Vector3d::Zero());
	for (auto face_index : current_node_mesh_rotated.faces()) {
		const std::size_t normal_index = static_cast<std::size_t>(face_index.idx());
		if (normal_index >= face_normals.size()) {
			continue;
		}
		Eigen::Vector3d normal;
		if (ComputeSelfSupportFaceNormal(
			current_node_mesh_rotated,
			face_index,
			normal)) {
			face_normals[normal_index] = normal;
		}
	}


	std::vector<SurfaceMesh::Vertex_index> vi_z_sorted;
	vi_z_sorted.reserve(current_node_mesh_rotated.number_of_vertices());
	for (auto v : current_node_mesh_rotated.vertices()) {
		vi_z_sorted.push_back(v);
	}
	std::sort(vi_z_sorted.begin(), vi_z_sorted.end(), [&](const SurfaceMesh::Vertex_index& a, const SurfaceMesh::Vertex_index& b) {
		return current_node_mesh_rotated.point(a).z() < current_node_mesh_rotated.point(b).z();
		});

	std::map<SurfaceMesh::Face_index, int> activeTriangles;

	double min_z = current_node_mesh_rotated.point(vi_z_sorted.front()).z();
	double layer_z = min_z + layer_height_value;
	int layer_height_adjust_count = 0;
	for (auto vi : vi_z_sorted) {
		double vi_z = current_node_mesh_rotated.point(vi).z();
		const auto& point = current_node_mesh_rotated.point(vi);
		//std::cout << "vi: " << vi << " z: " << vi_z << std::endl;

		while (layer_z <= vi_z) {
			if (vi_z - layer_z < layer_vertex_threshold) {
				layer_z += small_layer_height;
				++layer_height_adjust_count;
				std::cout << "Adjusting layer height for vertex at z = " << vi_z << " layer_height_adjust_count: " << layer_height_adjust_count << std::endl;
				continue;
			}


			SurfaceMeshSliceData slice_data;
			if (!BuildSurfaceMeshSingleSlice(activeTriangles, slice_data, layer_z)) {
				layer_z += small_layer_height;
				++layer_height_adjust_count;
				std::cout << "Adjusting layer height for vertex at z = " << vi_z << " layer_height_adjust_count: " << layer_height_adjust_count << std::endl;
				continue;
			}

			// build slice successfully
			layer_height_adjust_count = 0;
			if (slice_data.contour_points.empty()) {
				std::cout << "empty slice at z = " << layer_z << std::endl;
			}
			else {
				// GetTrianglesForSurfaceMeshSlices concatenates these arrays before
				// indexing them with global face IDs, so one complete copy is enough.
				if (slices.empty()) {
					slice_data.face_normals = face_normals;
				}
				slices.push_back(std::move(slice_data));
			}
			layer_z += layer_height_value;
		}

		// Update activeTriangles
		auto hf_vi = current_node_mesh_rotated.halfedge(vi);
		for (auto hf : current_node_mesh_rotated.halfedges_around_target(hf_vi)) {
			auto face_index = current_node_mesh_rotated.face(hf);
			if (face_index == SurfaceMesh::null_face()) {
				std::cout << "[HybridManufacturing::BuildSurfaceMeshSlices] Vertex " << vi << " has a halfedge with null face." << std::endl;
				continue;
			}
			auto active_it = activeTriangles.find(face_index);
			if (active_it == activeTriangles.end()) {
				activeTriangles[face_index] = 1;
			}
			else {
				active_it->second += 1;
				if (active_it->second == 3) {
					activeTriangles.erase(active_it);
				}
			}
		}

	}
	if (activeTriangles.size() > 0) {
		std::cout << "[HybridManufacturing::BuildSurfaceMeshSlices] Warning: " << activeTriangles.size() << " active triangles remain after processing all vertices." << std::endl;
	}
}

bool HybridManufacturing::BuildSurfaceMeshSingleSlice(std::map<SurfaceMesh::Face_index, int>& activeTriangles, SurfaceMeshSliceData& slice_data, double layer_z) const
{
	if (activeTriangles.empty()) {
		return true;
	}

	std::vector<bool> face_visited(activeTriangles.size(), false);
	std::vector<SurfaceMesh::Face_index> active_face_indices;
	active_face_indices.reserve(activeTriangles.size());
	std::unordered_map<SurfaceMesh::Face_index, int> face_index_to_active_index;
	for (auto active_triangle : activeTriangles) {
		const auto& face_index = active_triangle.first;
		active_face_indices.push_back(face_index);
		face_index_to_active_index[face_index] = active_face_indices.size() - 1;
	}

	for (size_t i = 0; i < active_face_indices.size(); ++i) {
		if (face_visited[i]) {
			continue;
		}
		const auto start_face_index = active_face_indices[i];
		if (start_face_index == SurfaceMesh::null_face()) {
			std::cout << "[HybridManufacturing::BuildSurfaceMeshSingleSlice] Encountered null face in activeTriangles." << std::endl;
			continue;
		}

		// Check if the face has vertices close to the slicing plane
		auto hf = current_node_mesh_rotated.halfedge(start_face_index);
		if (hf == SurfaceMesh::null_halfedge()) {
			std::cout << "[HybridManufacturing::BuildSurfaceMeshSingleSlice] Face " << start_face_index << " has null halfedge." << std::endl;
			continue;
		}


		for (int j = 0; j < 3; j++) {
			auto v = current_node_mesh_rotated.target(hf);
			const auto& p = current_node_mesh_rotated.point(v);
			if (std::abs(p.z() - layer_z) < layer_vertex_threshold) {
				return false;
			}
			hf = current_node_mesh_rotated.next(hf);
		}
		if (hf != current_node_mesh_rotated.halfedge(start_face_index)) {
			std::cout << "[HybridManufacturing::BuildSurfaceMeshSingleSlice] Face " << start_face_index << " does not have 3 vertices." << std::endl;
			return false;
		}

		// Find the intersection of the face with the slicing plane
		bool intersected = false;
		for (int j = 0; j < 3; j++) {
			auto v0 = current_node_mesh_rotated.source(hf);
			auto v1 = current_node_mesh_rotated.target(hf);
			const auto& p0 = current_node_mesh_rotated.point(v0);
			const auto& p1 = current_node_mesh_rotated.point(v1);
			if (!(p0.z() < layer_z && p1.z() > layer_z)) {
				hf = current_node_mesh_rotated.next(hf);
				continue;
			}
			intersected = true;
			Point_3 intersection_point = p0 + (layer_z - p0.z()) / (p1.z() - p0.z()) * (p1 - p0);
			slice_data.contour_points.push_back(Polyline_type{ intersection_point });
			slice_data.contour_face_ids.push_back({ start_face_index });
			slice_data.layer_z = layer_z;

			break;
		}

		if (!intersected) {
			std::cout << "[HybridManufacturing::BuildSurfaceMeshSingleSlice] Face " << start_face_index << " does not intersect with the slicing plane. Vertices: ";
			if (start_face_index == SurfaceMesh::null_face()) {
				std::cout << "null face" << std::endl;
				return false;
			}
			auto debug_hf = current_node_mesh_rotated.halfedge(start_face_index);
			for (int k = 0; k < 3; ++k) {
				const auto& p = current_node_mesh_rotated.point(current_node_mesh_rotated.target(debug_hf));
				std::cout << "(" << p.x() << ", " << p.y() << ", " << p.z() << ")";
				if (k < 2) {
					std::cout << ", ";
				}
				debug_hf = current_node_mesh_rotated.next(debug_hf);
			}
			std::cout << std::endl;
			std::cout << "layer_z: " << layer_z << std::endl;
			return false;
		}

		face_visited[i] = true;

		const auto start_halfedge = hf;

		hf = current_node_mesh_rotated.opposite(hf);
		while (current_node_mesh_rotated.face(hf) != start_face_index) {
			// Find the intersection of the face with the slicing plane
			bool intersected = false;
			for (int j = 0; j < 2; j++) {
				hf = current_node_mesh_rotated.next(hf);
				auto v0 = current_node_mesh_rotated.source(hf);
				auto v1 = current_node_mesh_rotated.target(hf);
				const auto& p0 = current_node_mesh_rotated.point(v0);
				const auto& p1 = current_node_mesh_rotated.point(v1);
				if (!(p0.z() < layer_z && p1.z() > layer_z)) {
					continue;
				}
				intersected = true;
				Point_3 intersection_point = p0 + (layer_z - p0.z()) / (p1.z() - p0.z()) * (p1 - p0);
				auto face_index = current_node_mesh_rotated.face(hf);
				auto active_index_it = face_index_to_active_index.find(face_index);

				if (active_index_it == face_index_to_active_index.end()) {
					std::cout << "[HybridManufacturing::BuildSurfaceMeshSingleSlice] Face " << face_index << " is not in activeTriangles." << std::endl;
					return false;
				}
				if (face_visited[active_index_it->second]) {
					std::cout << "[HybridManufacturing::BuildSurfaceMeshSingleSlice] Face " << face_index << " has already been visited." << std::endl;
					return false;
				}
				face_visited[active_index_it->second] = true;


				slice_data.contour_points.back().push_back(intersection_point);
				slice_data.contour_face_ids.back().push_back(face_index);

				break;
			}

			if (!intersected) {
				auto face_index = current_node_mesh_rotated.face(hf);
				std::cout << "[HybridManufacturing::BuildSurfaceMeshSingleSlice] Face " << face_index << " does not intersect with the slicing plane." << std::endl;
				if (face_index == SurfaceMesh::null_face()) {
					std::cout << "[HybridManufacturing::BuildSurfaceMeshSingleSlice] Halfedge " << hf << " has null face." << std::endl;
					return false;
				}
				auto debug_hf = current_node_mesh_rotated.halfedge(face_index);
				if (debug_hf == SurfaceMesh::null_halfedge()) {
					std::cout << "[HybridManufacturing::BuildSurfaceMeshSingleSlice] Face " << face_index << " has null halfedge." << std::endl;
					return false;
				}
				for (int k = 0; k < 3; ++k) {
					const auto& p = current_node_mesh_rotated.point(current_node_mesh_rotated.target(debug_hf));
					std::cout << "(" << p.x() << ", " << p.y() << ", " << p.z() << ")";
					if (k < 2) {
						std::cout << ", ";
					}
					debug_hf = current_node_mesh_rotated.next(debug_hf);
				}
				std::cout << std::endl;
				std::cout << "layer_z: " << layer_z << std::endl;
				return false;
				return false;
			}

			hf = current_node_mesh_rotated.opposite(hf);
		}
	}

	slice_data.contour_parent_ids =
		vasco::contour_containment::BuildContourParentIds(slice_data.contour_points);
	return true;
}

namespace
{
	using CutEdgeKey = std::pair<int, int>;

	CutEdgeKey MakeCutEdgeKey(int a, int b)
	{
		if (a > b) {
			std::swap(a, b);
		}
		return { a, b };
	}

	std::uint64_t EncodeCutEdge(int a, int b)
	{
		if (a > b) {
			std::swap(a, b);
		}
		return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(a)) << 32)
			| static_cast<std::uint32_t>(b);
	}

	bool TryGetUpperCutBoundaryEdge(
		const Slicer_2& slicer,
		const TRiangle& triangle,
		double plane_z,
		CutEdgeKey& boundary_edge,
		double eps = 1e-4)
	{
		std::array<int, 2> plane_vertices{};
		int plane_vertex_count = 0;
		int off_plane_vertex = -1;
		for (int k = 0; k < 3; ++k) {
			const int vertex_id = triangle[k];
			if (vertex_id < 0 || vertex_id >= static_cast<int>(slicer.positions.size())) {
				return false;
			}
			const double z = slicer.positions[vertex_id][2];
			if (std::abs(z - plane_z) <= eps) {
				if (plane_vertex_count < 2) {
					plane_vertices[plane_vertex_count] = vertex_id;
				}
				++plane_vertex_count;
			}
			else {
				off_plane_vertex = vertex_id;
			}
		}

		if (plane_vertex_count != 2 || off_plane_vertex < 0) {
			return false;
		}
		if (slicer.positions[off_plane_vertex][2] <= plane_z + eps) {
			return false;
		}

		boundary_edge = MakeCutEdgeKey(plane_vertices[0], plane_vertices[1]);
		return true;
	}

	bool CutBoundaryEdgeBelongsToLayer(
		const Slicer_2& slicer,
		const CutEdgeKey& edge,
		const LayerContainmentPolygon& layer_polygon,
		double tolerance)
	{
		const auto& a = slicer.positions[edge.first];
		const auto& b = slicer.positions[edge.second];
		const Point_2 midpoint((a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5);
		return IsPointInsidePreparedLayerPolygon(layer_polygon, midpoint, tolerance);
	}

	double SquaredDistanceCutEdgeToLayerBoundary(
		const Slicer_2& slicer,
		const CutEdgeKey& edge,
		const LayerContainmentPolygon& layer_polygon)
	{
		if (layer_polygon.points.size() < 2) {
			return std::numeric_limits<double>::max();
		}
		const auto& a = slicer.positions[edge.first];
		const auto& b = slicer.positions[edge.second];
		const Point_2 midpoint((a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5);
		double min_distance_sq = std::numeric_limits<double>::max();
		for (std::size_t i = 0; i < layer_polygon.points.size(); ++i) {
			min_distance_sq = std::min(
				min_distance_sq,
				SquaredDistancePointToSegment2D(
					midpoint,
					layer_polygon.points[i],
					layer_polygon.points[(i + 1) % layer_polygon.points.size()]));
		}
		return min_distance_sq;
	}

	std::vector<int> CollectCutRemovalSeedFaces(
		const Slicer_2& slicer,
		const std::vector<std::vector<cv::Point3d>>& cut_layers,
		std::vector<std::vector<std::pair<int, int>>>& cutting_plane_edges)
	{
		cutting_plane_edges.assign(cut_layers.size(), {});
		std::vector<LayerContainmentPolygon> layer_polygons;
		layer_polygons.reserve(cut_layers.size());
		for (const auto& layer : cut_layers) {
			std::vector<Point_2> points;
			points.reserve(layer.size());
			for (const auto& point : layer) {
				points.emplace_back(point.x, point.y);
			}
			layer_polygons.push_back(PrepareLayerContainmentPolygon(points, 1e-8));
		}

		std::vector<bool> is_seed_face(slicer.triangles.size(), false);
		std::vector<int> seed_face_ids;
		std::vector<std::unordered_set<std::uint64_t>> layer_edge_keys(cut_layers.size());
		std::vector<std::pair<double, int>> layers_by_height;
		layers_by_height.reserve(cut_layers.size());
		for (int layer_id = 0; layer_id < static_cast<int>(cut_layers.size()); ++layer_id) {
			if (!cut_layers[layer_id].empty() && layer_polygons[layer_id].points.size() >= 3) {
				layers_by_height.emplace_back(cut_layers[layer_id][0].z, layer_id);
			}
		}
		std::sort(layers_by_height.begin(), layers_by_height.end());

		constexpr double plane_eps = 1e-4;
		for (int face_id = 0; face_id < static_cast<int>(slicer.triangles.size()); ++face_id) {
			const auto& triangle = slicer.triangles[face_id];
			for (int edge_index = 0; edge_index < 3; ++edge_index) {
				const int a_id = triangle[edge_index];
				const int b_id = triangle[(edge_index + 1) % 3];
				const int opposite_id = triangle[(edge_index + 2) % 3];
				const auto& a = slicer.positions[a_id];
				const auto& b = slicer.positions[b_id];
				const auto& opposite = slicer.positions[opposite_id];
				if (std::abs(a[2] - b[2]) > plane_eps) {
					continue;
				}
				const double edge_z = (a[2] + b[2]) * 0.5;
				if (opposite[2] <= edge_z + plane_eps) {
					continue;
				}

				const auto first_layer = std::lower_bound(
					layers_by_height.begin(),
					layers_by_height.end(),
					std::make_pair(edge_z - plane_eps, std::numeric_limits<int>::min()));
				int closest_layer_id = -1;
				double closest_layer_distance_sq = std::numeric_limits<double>::max();
				const CutEdgeKey boundary_edge = MakeCutEdgeKey(a_id, b_id);
				for (auto layer_it = first_layer;
					layer_it != layers_by_height.end() && layer_it->first <= edge_z + plane_eps;
					++layer_it) {
					const int layer_id = layer_it->second;
					const double distance_sq = SquaredDistanceCutEdgeToLayerBoundary(
						slicer,
						boundary_edge,
						layer_polygons[layer_id]);
					if (distance_sq < closest_layer_distance_sq) {
						closest_layer_distance_sq = distance_sq;
						closest_layer_id = layer_id;
					}
				}

				if (closest_layer_id < 0 || closest_layer_distance_sq > 1e-4) {
					continue;
				}
				const std::uint64_t edge_key = EncodeCutEdge(boundary_edge.first, boundary_edge.second);
				if (layer_edge_keys[closest_layer_id].insert(edge_key).second) {
					cutting_plane_edges[closest_layer_id].push_back(boundary_edge);
				}
				if (!is_seed_face[face_id]) {
					is_seed_face[face_id] = true;
					seed_face_ids.push_back(face_id);
				}
			}
		}
		return seed_face_ids;
	}

	std::vector<bool> ExpandRemovedFacesAcrossNonCutEdges(
		const Slicer_2& slicer,
		const std::vector<int>& seed_face_ids,
		const std::vector<std::vector<std::pair<int, int>>>& cutting_plane_edges)
	{
		std::unordered_set<std::uint64_t> barrier_edges;
		std::size_t barrier_edge_count = 0;
		for (const auto& layer_edges : cutting_plane_edges) {
			barrier_edge_count += layer_edges.size();
		}
		barrier_edges.reserve(barrier_edge_count * 2 + 1);
		for (const auto& layer_edges : cutting_plane_edges) {
			for (const auto& edge : layer_edges) {
				barrier_edges.insert(EncodeCutEdge(edge.first, edge.second));
			}
		}

		std::unordered_map<std::uint64_t, std::vector<int>> edge_to_faces;
		edge_to_faces.reserve(slicer.triangles.size() * 2);
		for (int face_id = 0; face_id < static_cast<int>(slicer.triangles.size()); ++face_id) {
			const auto& tri = slicer.triangles[face_id];
			for (int k = 0; k < 3; ++k) {
				edge_to_faces[EncodeCutEdge(tri[k], tri[(k + 1) % 3])].push_back(face_id);
			}
		}

		std::vector<bool> selected(slicer.triangles.size(), false);
		std::queue<int> frontier;
		for (int face_id : seed_face_ids) {
			if (face_id < 0 || face_id >= static_cast<int>(selected.size()) || selected[face_id]) {
				continue;
			}
			selected[face_id] = true;
			frontier.push(face_id);
		}

		while (!frontier.empty()) {
			const int face_id = frontier.front();
			frontier.pop();
			const auto& tri = slicer.triangles[face_id];
			for (int k = 0; k < 3; ++k) {
				const std::uint64_t edge = EncodeCutEdge(tri[k], tri[(k + 1) % 3]);
				if (barrier_edges.find(edge) != barrier_edges.end()) {
					continue;
				}
				const auto face_it = edge_to_faces.find(edge);
				if (face_it == edge_to_faces.end()) {
					continue;
				}
				for (int neighbor_face_id : face_it->second) {
					if (!selected[neighbor_face_id]) {
						selected[neighbor_face_id] = true;
						frontier.push(neighbor_face_id);
					}
				}
			}
		}

		return selected;
	}

	// Finds unresolved subtractive-inaccessible source vertices referenced by
	// the final removed surface. Quantized lookup keeps this validation near
	// linearithmic instead of comparing every area_S with every removed face.
	std::vector<std::pair<int, int>> FindUnresolvedAreaSOnRemovedSurface(
		const Slicer_2& slicer,
		const std::vector<TRiangle>& removed_surface_triangles,
		const Eigen::MatrixXd& original_vertices,
		const std::unordered_map<int, int>& map_s_to_vertex,
		const std::vector<bool>& searched_s_flags)
	{
		std::set<vasco::contact_triangulation::QuantizedPointKey> removed_vertex_keys;
		for (const TRiangle& triangle : removed_surface_triangles) {
			for (int corner = 0; corner < 3; ++corner) {
				const int vertex_id = triangle[corner];
				if (vertex_id < 0 || vertex_id >= static_cast<int>(slicer.positions.size())) {
					continue;
				}
				const auto& point = slicer.positions[vertex_id];
				removed_vertex_keys.insert(vasco::contact_triangulation::MakeQuantizedKey(
					Point_3(point[0], point[1], point[2]),
					1e6));
			}
		}

		std::vector<std::pair<int, int>> hits;
		for (const auto& entry : map_s_to_vertex) {
			const int s_id = entry.first;
			const int original_vertex_id = entry.second;
			if (s_id < 0
				|| s_id >= static_cast<int>(searched_s_flags.size())
				|| searched_s_flags[s_id]
				|| original_vertex_id < 0
				|| original_vertex_id >= original_vertices.rows()) {
				continue;
			}

			const auto key = vasco::contact_triangulation::MakeQuantizedKey(
				Point_3(
					original_vertices(original_vertex_id, 0),
					original_vertices(original_vertex_id, 1),
					original_vertices(original_vertex_id, 2)),
				1e6);
			if (removed_vertex_keys.find(key) != removed_vertex_keys.end()) {
				hits.emplace_back(s_id, original_vertex_id);
			}
		}
		std::sort(hits.begin(), hits.end());
		return hits;
	}

	std::vector<std::vector<int>> CollectUnselectedFaceComponents(
		const Slicer_2& slicer,
		const std::vector<bool>& selected_faces)
	{
		const int face_count = static_cast<int>(slicer.triangles.size());
		std::vector<int> parent(face_count, -1);
		std::vector<int> rank(face_count, 0);
		for (int face_id = 0; face_id < face_count; ++face_id) {
			if (face_id >= static_cast<int>(selected_faces.size()) || !selected_faces[face_id]) {
				parent[face_id] = face_id;
			}
		}

		auto find_root = [&](int face_id) {
			int root = face_id;
			while (parent[root] != root) {
				root = parent[root];
			}
			while (parent[face_id] != face_id) {
				const int next = parent[face_id];
				parent[face_id] = root;
				face_id = next;
			}
			return root;
		};

		auto unite = [&](int a, int b) {
			int root_a = find_root(a);
			int root_b = find_root(b);
			if (root_a == root_b) {
				return;
			}
			if (rank[root_a] < rank[root_b]) {
				std::swap(root_a, root_b);
			}
			parent[root_b] = root_a;
			if (rank[root_a] == rank[root_b]) {
				++rank[root_a];
			}
		};

		std::unordered_map<std::uint64_t, int> first_face_for_edge;
		first_face_for_edge.reserve(slicer.triangles.size() * 2);
		for (int face_id = 0; face_id < face_count; ++face_id) {
			if (parent[face_id] < 0) {
				continue;
			}
			const auto& tri = slicer.triangles[face_id];
			for (int k = 0; k < 3; ++k) {
				const std::uint64_t edge = EncodeCutEdge(tri[k], tri[(k + 1) % 3]);
				const auto [it, inserted] = first_face_for_edge.emplace(edge, face_id);
				if (!inserted) {
					unite(face_id, it->second);
				}
			}
		}

		std::unordered_map<int, std::vector<int>> faces_by_root;
		faces_by_root.reserve(face_count / 8 + 1);
		for (int face_id = 0; face_id < face_count; ++face_id) {
			if (parent[face_id] >= 0) {
				faces_by_root[find_root(face_id)].push_back(face_id);
			}
		}

		std::vector<std::vector<int>> components;
		components.reserve(faces_by_root.size());
		for (auto& [root, faces] : faces_by_root) {
			components.push_back(std::move(faces));
		}
		return components;
	}

	std::vector<std::vector<std::vector<int>>> ExtractAllClosedCutBoundaryLoops(
		const std::vector<std::vector<std::pair<int, int>>>& cutting_plane_edges,
		bool& all_loops_valid)
	{
		all_loops_valid = true;
		std::vector<std::vector<std::vector<int>>> loops_by_layer(cutting_plane_edges.size());
		for (int layer_id = 0; layer_id < static_cast<int>(cutting_plane_edges.size()); ++layer_id) {
			std::set<CutEdgeKey> unused_edges;
			std::map<int, std::vector<int>> vertex_neighbors;
			for (const auto& raw_edge : cutting_plane_edges[layer_id]) {
				const CutEdgeKey edge = MakeCutEdgeKey(raw_edge.first, raw_edge.second);
				if (edge.first == edge.second || !unused_edges.insert(edge).second) {
					continue;
				}
				vertex_neighbors[edge.first].push_back(edge.second);
				vertex_neighbors[edge.second].push_back(edge.first);
			}

			bool layer_valid = !unused_edges.empty();
			for (const auto& [vertex_id, neighbors] : vertex_neighbors) {
				if (neighbors.size() != 2) {
					std::cout << "[HybridManufacturing::CutMesh] Cut boundary vertex degree is not 2: "
						<< "layer=" << layer_id
						<< ", vertex=" << vertex_id
						<< ", degree=" << neighbors.size()
						<< std::endl;
					layer_valid = false;
				}
			}

			if (!layer_valid) {
				all_loops_valid = false;
				continue;
			}

			while (!unused_edges.empty()) {
				const int start_vertex = unused_edges.begin()->first;
				int current_vertex = start_vertex;
				std::vector<int> loop;
				bool closed = false;
				for (std::size_t step = 0; step <= vertex_neighbors.size(); ++step) {
					loop.push_back(current_vertex);
					int next_vertex = -1;
					for (int neighbor : vertex_neighbors[current_vertex]) {
						const CutEdgeKey edge = MakeCutEdgeKey(current_vertex, neighbor);
						if (unused_edges.find(edge) != unused_edges.end()) {
							next_vertex = neighbor;
							unused_edges.erase(edge);
							break;
						}
					}
					if (next_vertex < 0) {
						break;
					}
					current_vertex = next_vertex;
					if (current_vertex == start_vertex) {
						closed = true;
						break;
					}
				}

				if (!closed || loop.size() < 3) {
					std::cout << "[HybridManufacturing::CutMesh] Failed to extract closed cut loop: layer="
						<< layer_id << std::endl;
					all_loops_valid = false;
					loops_by_layer[layer_id].clear();
					break;
				}
				loops_by_layer[layer_id].push_back(std::move(loop));
			}
		}
		return loops_by_layer;
	}

	struct RemainingCutBoundaryLoopStats {
		std::size_t input_loop_count = 0;
		std::size_t retained_loop_count = 0;
		std::size_t discarded_orphan_loop_count = 0;
		std::size_t partial_boundary_loop_count = 0;
	};

	// Keeps only cut loops that are actual open boundaries of the remaining mesh.
	// A loop with no remaining boundary edges belongs entirely to the removed side and
	// must not be capped; a partially matching loop indicates inconsistent face removal.
	std::vector<std::vector<std::vector<int>>> FilterCutLoopsToRemainingMeshBoundary(
		const Slicer_2& remaining_mesh,
		const std::vector<std::vector<std::vector<int>>>& loops_by_layer,
		RemainingCutBoundaryLoopStats& stats,
		bool& all_loops_valid)
	{
		stats = {};
		std::unordered_map<std::uint64_t, std::size_t> remaining_edge_face_counts;
		remaining_edge_face_counts.reserve(remaining_mesh.triangles.size() * 2 + 1);
		for (const auto& triangle : remaining_mesh.triangles) {
			for (int edge_index = 0; edge_index < 3; ++edge_index) {
				++remaining_edge_face_counts[EncodeCutEdge(
					triangle[edge_index],
					triangle[(edge_index + 1) % 3])];
			}
		}

		std::vector<std::vector<std::vector<int>>> filtered_loops(loops_by_layer.size());
		for (std::size_t layer_id = 0; layer_id < loops_by_layer.size(); ++layer_id) {
			for (const auto& loop : loops_by_layer[layer_id]) {
				++stats.input_loop_count;
				std::size_t remaining_boundary_edge_count = 0;
				for (std::size_t vertex_index = 0; vertex_index < loop.size(); ++vertex_index) {
					const std::uint64_t edge = EncodeCutEdge(
						loop[vertex_index],
						loop[(vertex_index + 1) % loop.size()]);
					const auto edge_it = remaining_edge_face_counts.find(edge);
					if (edge_it != remaining_edge_face_counts.end() && edge_it->second == 1) {
						++remaining_boundary_edge_count;
					}
				}

				if (remaining_boundary_edge_count == loop.size()) {
					filtered_loops[layer_id].push_back(loop);
					++stats.retained_loop_count;
					continue;
				}
				if (remaining_boundary_edge_count == 0) {
					++stats.discarded_orphan_loop_count;
					continue;
				}

				++stats.partial_boundary_loop_count;
				all_loops_valid = false;
				std::cerr << "[HybridManufacturing::CutMesh] Cut loop only partially belongs to "
					<< "the remaining mesh boundary: layer=" << layer_id
					<< ", boundary_edges=" << remaining_boundary_edge_count
					<< ", loop_edges=" << loop.size()
					<< std::endl;
			}
		}
		return filtered_loops;
	}

	double CutRingSignedArea(const std::vector<int>& ring, const Slicer_2& slicer)
	{
		double area_twice = 0.0;
		for (std::size_t i = 0; i < ring.size(); ++i) {
			const auto& a = slicer.positions[ring[i]];
			const auto& b = slicer.positions[ring[(i + 1) % ring.size()]];
			area_twice += a[0] * b[1] - b[0] * a[1];
		}
		return area_twice * 0.5;
	}

	bool CutRingContainsPoint(
		const std::vector<int>& ring,
		const Slicer_2& slicer,
		const Slicer_2::Vec3& point)
	{
		bool inside = false;
		for (std::size_t i = 0, j = ring.size() - 1; i < ring.size(); j = i++) {
			const auto& a = slicer.positions[ring[i]];
			const auto& b = slicer.positions[ring[j]];
			const bool crosses = ((a[1] > point[1]) != (b[1] > point[1]))
				&& (point[0] < (b[0] - a[0]) * (point[1] - a[1]) / (b[1] - a[1] + 1e-30) + a[0]);
			if (crosses) {
				inside = !inside;
			}
		}
		return inside;
	}

	inline double PointDistance2D(const Point_3& a, const Point_3& b)
	{
		double dx = a.x() - b.x();
		double dy = a.y() - b.y();
		return std::sqrt(dx * dx + dy * dy);
	}

	inline bool PointsClose(const Point_3& a, const Point_3& b, double tol)
	{
		return PointDistance2D(a, b) <= tol;
	}

	struct ReferencedSlicerZBounds {
		double min_z = 0.0;
		double max_z = 0.0;
		std::size_t referenced_vertex_count = 0;
		std::size_t invalid_vertex_index_count = 0;
	};

	bool ComputeReferencedSlicerZBounds(
		const Slicer_2& slicer,
		ReferencedSlicerZBounds& bounds)
	{
		bounds = {};
		std::vector<bool> referenced_vertices(slicer.positions.size(), false);
		bool has_finite_referenced_vertex = false;

		for (const auto& triangle : slicer.triangles) {
			for (int corner = 0; corner < 3; ++corner) {
				const int vertex_id = triangle[corner];
				if (vertex_id < 0 || vertex_id >= static_cast<int>(slicer.positions.size())) {
					++bounds.invalid_vertex_index_count;
					continue;
				}
				if (referenced_vertices[vertex_id]) {
					continue;
				}

				referenced_vertices[vertex_id] = true;
				++bounds.referenced_vertex_count;
				const double z = slicer.positions[vertex_id][2];
				if (!std::isfinite(z)) {
					continue;
				}

				if (!has_finite_referenced_vertex) {
					bounds.min_z = z;
					bounds.max_z = z;
					has_finite_referenced_vertex = true;
				}
				else {
					bounds.min_z = std::min(bounds.min_z, z);
					bounds.max_z = std::max(bounds.max_z, z);
				}
			}
		}

		return has_finite_referenced_vertex;
	}

	struct OriginalFacePlane {
		Eigen::Vector3d origin = Eigen::Vector3d::Zero();
		Eigen::Vector3d unit_normal = Eigen::Vector3d::Zero();
		bool is_valid = false;
	};

	// Computes a model-scale tolerance for deciding whether a candidate face
	// actually lies on an original input face instead of merely projecting onto it.
	double ComputeOriginalSurfacePlaneTolerance(const Eigen::MatrixXd& vertices)
	{
		if (vertices.rows() == 0 || vertices.cols() < 3) {
			return 1e-6;
		}

		const Eigen::RowVector3d minimum = vertices.leftCols<3>().colwise().minCoeff();
		const Eigen::RowVector3d maximum = vertices.leftCols<3>().colwise().maxCoeff();
		const double model_diagonal = (maximum - minimum).norm();
		if (!std::isfinite(model_diagonal)) {
			return 1e-6;
		}
		return std::clamp(model_diagonal * 1e-6, 1e-6, 1e-3);
	}

	// Precomputes original face planes once so the surface filter can reject
	// historical cap faces without repeatedly rebuilding plane normals.
	std::vector<OriginalFacePlane> BuildOriginalFacePlanes(
		const Eigen::MatrixXd& vertices,
		const Eigen::MatrixXi& faces,
		std::size_t face_count)
	{
		std::vector<OriginalFacePlane> planes(face_count);
		for (std::size_t face_id = 0; face_id < face_count; ++face_id) {
			const int v0 = faces(static_cast<int>(face_id), 0);
			const int v1 = faces(static_cast<int>(face_id), 1);
			const int v2 = faces(static_cast<int>(face_id), 2);
			if (v0 < 0 || v1 < 0 || v2 < 0
				|| v0 >= vertices.rows()
				|| v1 >= vertices.rows()
				|| v2 >= vertices.rows()) {
				continue;
			}

			OriginalFacePlane& plane = planes[face_id];
			plane.origin = vertices.row(v0).transpose();
			const Eigen::Vector3d edge01 =
				(vertices.row(v1) - vertices.row(v0)).transpose();
			const Eigen::Vector3d edge02 =
				(vertices.row(v2) - vertices.row(v0)).transpose();
			plane.unit_normal = edge01.cross(edge02);
			const double normal_length = plane.unit_normal.norm();
			if (!plane.unit_normal.allFinite() || normal_length <= 1e-12) {
				plane.unit_normal.setZero();
				continue;
			}
			plane.unit_normal /= normal_length;
			plane.is_valid = true;
		}
		return planes;
	}

	// Returns true only when every candidate vertex is within the allowed
	// perpendicular distance of the original face plane.
	bool IsTriangleOnOriginalFacePlane(
		const Eigen::Vector3d& v0,
		const Eigen::Vector3d& v1,
		const Eigen::Vector3d& v2,
		const OriginalFacePlane& plane,
		double tolerance)
	{
		if (!plane.is_valid) {
			return false;
		}
		return std::abs((v0 - plane.origin).dot(plane.unit_normal)) <= tolerance
			&& std::abs((v1 - plane.origin).dot(plane.unit_normal)) <= tolerance
			&& std::abs((v2 - plane.origin).dot(plane.unit_normal)) <= tolerance;
	}
}

HybridManufacturing::HybridManufacturing(std::string input_folder, std::string suf, Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::MatrixXd N)
{
	this->V = V;
	this->F = F;
	this->Normals = N;
	this->file_name = input_folder + "\\" + suf + "\\" + suf;
	this->input_folder = input_folder + "\\" + suf;
	this->suf = suf;
	open_vis_additive_accessibility_debug = false;
	open_vis_additive_self_support_debug = false;
	InitializeSurfaceMeshFromVF();
	InitializePolyscope();
}

HybridManufacturing::~HybridManufacturing()
{
}

std::string HybridManufacturing::VisDir() const
{
	const std::string dir = input_folder + "\\vis";
	EnsureDirectory(dir);
	return dir;
}

std::string HybridManufacturing::VisPath(const std::string& output_file_name) const
{
	return JoinPath(VisDir(), output_file_name);
}

// Returns the per-model directory for accessibility and inaccessible-point diagnostics.
std::string HybridManufacturing::UnaccessibleVisDir() const
{
	const std::string dir = JoinPath(VisDir(), "unaccessible");
	EnsureDirectory(dir);
	return dir;
}

// Builds a path under vis\unaccessible and creates the category directory on demand.
std::string HybridManufacturing::UnaccessibleVisPath(const std::string& output_file_name) const
{
	return JoinPath(UnaccessibleVisDir(), output_file_name);
}

// Returns the per-model directory for subtractive patch outputs and merged-patch diagnostics.
std::string HybridManufacturing::PatchVisDir() const
{
	const std::string dir = JoinPath(VisDir(), "patch");
	EnsureDirectory(dir);
	return dir;
}

// Builds a path under vis\patch and creates the category directory on demand.
std::string HybridManufacturing::PatchVisPath(const std::string& output_file_name) const
{
	return JoinPath(PatchVisDir(), output_file_name);
}

// Builds a path under the dedicated ancestor-source report directory.
std::string HybridManufacturing::AncestorSourceVisPath(const std::string& output_file_name) const
{
	const std::string dir = JoinPath(VisDir(), "ancestor_sources");
	EnsureDirectory(dir);
	return JoinPath(dir, output_file_name);
}

void HybridManufacturing::SetContactFaceTriangulationMode(ContactFaceTriangulationMode mode)
{
	contact_face_triangulation_mode = mode;
}

HybridManufacturing::ContactFaceTriangulationMode HybridManufacturing::GetContactFaceTriangulationMode() const
{
	return contact_face_triangulation_mode;
}

void HybridManufacturing::InitializeSurfaceMeshFromVF()
{
	input_mesh.clear();
	if (V.rows() == 0 || F.rows() == 0) {
		return;
	}

	std::vector<SurfaceMesh::Vertex_index> vertices;
	vertices.reserve(V.rows());
	for (int i = 0; i < V.rows(); ++i) {
		vertices.push_back(input_mesh.add_vertex(Point_3(V(i, 0), V(i, 1), V(i, 2))));
	}

	for (int i = 0; i < F.rows(); ++i) {
		input_mesh.add_face(vertices[F(i, 0)], vertices[F(i, 1)], vertices[F(i, 2)]);
	}
}

std::vector<TRiangle> HybridManufacturing::FilterSurfaceRemoveTriangles(
	const Slicer_2& slicer,
	const std::vector<TRiangle>& remove_triangles) const
{
	std::vector<TRiangle> surface_triangles;
	surface_triangles.reserve(remove_triangles.size());

	const std::size_t face_cnt = static_cast<std::size_t>(
		std::min(F.rows(), Normals.rows()));
	const double bottom_filter_threshold = 2.0;
	const double bottom_normal_abs_z_threshold = 0.9;
	const double surface_plane_tolerance = ComputeOriginalSurfacePlaneTolerance(V);
	const std::vector<OriginalFacePlane> original_face_planes =
		BuildOriginalFacePlanes(V, F, face_cnt);
	double model_min_z = MAX_D;
	for (int i = 0; i < V.rows(); ++i) {
		model_min_z = std::min(model_min_z, V(i, 2));
	}

	std::size_t skipped_bottom_triangles = 0;
	for (size_t i = 0; i < remove_triangles.size(); i++) {
		Eigen::Vector3d v1, v2, v3;
		v1.x() = slicer.positions[remove_triangles[i][0]][0];
		v1.y() = slicer.positions[remove_triangles[i][0]][1];
		v1.z() = slicer.positions[remove_triangles[i][0]][2];
		v2.x() = slicer.positions[remove_triangles[i][1]][0];
		v2.y() = slicer.positions[remove_triangles[i][1]][1];
		v2.z() = slicer.positions[remove_triangles[i][1]][2];
		v3.x() = slicer.positions[remove_triangles[i][2]][0];
		v3.y() = slicer.positions[remove_triangles[i][2]][1];
		v3.z() = slicer.positions[remove_triangles[i][2]][2];
		double ans = (v2.x() - v1.x()) * (v2.y() - v3.y()) - (v2.y() - v1.y()) * (v2.x() - v3.x());
		if (ans > 0) {
			swap(v2, v3);
		}
		double na = (v2.y() - v1.y()) * (v3.z() - v1.z()) - (v2.z() - v1.z()) * (v3.y() - v1.y());
		double nb = (v2.z() - v1.z()) * (v3.x() - v1.x()) - (v2.x() - v1.x()) * (v3.z() - v1.z());
		double nc = (v2.x() - v1.x()) * (v3.y() - v1.y()) - (v2.y() - v1.y()) * (v3.x() - v1.x());

		Eigen::Vector3d normal_vector(na, nb, nc);
		const double normal_length = normal_vector.norm();
		if (normal_length <= 1e-12) {
			continue;
		}
		normal_vector /= normal_length;

		bool jud_surface = false;
		for (size_t j = 0; j < face_cnt; j++) {
			double n_dot = normal_vector.dot(Normals.row(j));
			if (fabs(n_dot) < 0.999) {
				continue;
			}
			if (!IsTriangleOnOriginalFacePlane(
				v1,
				v2,
				v3,
				original_face_planes[j],
				surface_plane_tolerance)) {
				continue;
			}
			Eigen::Vector3d v1_s(v1.x(), v1.y(), v1.z());
			Eigen::Vector3d v2_s(v2.x(), v2.y(), v2.z());
			Eigen::Vector3d v3_s(v3.x(), v3.y(), v3.z());

			bool b1 = vasco::core::isPointInTriangle(v1_s, V.row(F(j, 0)), V.row(F(j, 1)), V.row(F(j, 2)));
			bool b2 = vasco::core::isPointInTriangle(v2_s, V.row(F(j, 0)), V.row(F(j, 1)), V.row(F(j, 2)));
			bool b3 = vasco::core::isPointInTriangle(v3_s, V.row(F(j, 0)), V.row(F(j, 1)), V.row(F(j, 2)));

			if (b1 && b2 && b3) {
				jud_surface = true;
				break;
			}
		}

		if (jud_surface) {
			const double centroid_z = (v1.z() + v2.z() + v3.z()) / 3.0;
			const bool is_bottom_triangle =
				centroid_z - model_min_z <= bottom_filter_threshold
				&& fabs(normal_vector.z()) >= bottom_normal_abs_z_threshold;
			if (is_bottom_triangle) {
				++skipped_bottom_triangles;
				continue;
			}
			surface_triangles.push_back(remove_triangles[i]);
		}
	}

	if (skipped_bottom_triangles != 0) {
		std::cout << "[HybridManufacturing::FilterSurfaceRemoveTriangles] skipped bottom surface triangles: "
			<< skipped_bottom_triangles
			<< ", threshold_z=" << bottom_filter_threshold
			<< ", min_z=" << model_min_z
			<< std::endl;
	}

	return surface_triangles;
}

void HybridManufacturing::InitializeVoronoi()
{

	std::vector<vasco::VoronoiCell> voronoiCells;
	std::vector<Eigen::Vector3d> bottomVertices;

	std::vector<vasco::VoronoiCell> voronoiCells1;
	std::vector<Eigen::Vector3d> bottomVertices1;
	vasco::BuildVoronoiCells(V, F, 2.0, voronoiCells, bottomVertices,
		open_vis_voronoi, file_name);

	vasco::BuildVoronoiCells(input_mesh, 2.0, voronoiCells1, bottomVertices1,
		open_vis_voronoi, file_name);

	size_t cell_cnt = voronoiCells.size();
	for (size_t i = 0; i < cell_cnt; i++) {
		if (voronoiCells[i].is_available != voronoiCells1[i].is_available) {
			std::cerr << "Error: Voronoi cell availability mismatch at index " << i << std::endl;
		}
		if (voronoiCells[i].is_available) {
			size_t point_cnt = voronoiCells[i].all_points_in_polygon.size();
			if (point_cnt != voronoiCells1[i].all_points_in_polygon.size()) {
				std::cerr << "Error: Voronoi cell point count mismatch at index " << i << std::endl;
			}
			//for (size_t j = 0; j < point_cnt; j++) {
			//	if ((voronoiCells[i].all_points_in_polygon[j] - voronoiCells1[i].all_points_in_polygon[j]).norm() > 1e-6) {
			//		std::cerr << "Error: Voronoi cell point mismatch at index " << i << ", point " << j << std::endl;
			//	}
			//}
			if (voronoiCells[i].site != voronoiCells1[i].site) {
				std::cerr << "Error: Voronoi cell site mismatch at index " << i << std::endl;
			}
			if (voronoiCells[i].adjacent_cells.size() != voronoiCells1[i].adjacent_cells.size()) {
				std::cerr << "Error: Voronoi cell adjacent cell count mismatch at index " << i << std::endl;
			}
		}
	}
	all_voronoi_cells = voronoiCells;
	V_bottom = bottomVertices;
}

void HybridManufacturing::InitializePolyscope()
{
	polyscope::init();
	polyscope::view::setUpDir(polyscope::UpDir::ZUp);
}

int HybridManufacturing::CollisionDetectionForSubtractiveManufacturing(cutter the_nozzle)
{
	clock_t start_time, end_time;
	clock_t start_time_2, end_time_2;
	clock_t start_time_3, end_time_3;
	clock_t start_time_4, end_time_4;
	double sum_time = 0;
	int test_cont_num = 0;

	start_time = clock();
	std::cout << "Doing collision detection for subtractive manufacturing......" << endl;
	string file_name_2 = file_name;
	int cont_num = 0;
	sampling_subtractive.OrientationSamplePoints();	//sampling_subtractive生成球面采样点
	Eigen::Matrix3d rotMatrix;
	std::vector<bool> flag_accessible_points;	//true值代表对应的voronoi单元在某个采样方向下可达
	std::vector<bool> flag_covering_points;	//true值代表对应的voronoi单元在某个采样方向下被喷嘴覆盖
	//double nozzle_par = (the_nozzle.upper_surface_r - the_nozzle.lowwer_surface_r) / the_nozzle.nozzle_H_half;

	flag_accessible_points.resize(all_voronoi_cells.size());
	flag_covering_points.resize(all_voronoi_cells.size());
	flag_voronoi_has_been_printed.resize(all_voronoi_cells.size());
	for (int i = 0; i < all_voronoi_cells.size(); i++) {
		flag_accessible_points[i] = false;
		flag_voronoi_has_been_printed[i] = true;
		flag_covering_points[i] = false;
	}

	for (int ori = 0; ori < sampling_subtractive.sample_points.size(); ori++) {	//枚举所有采样方向

		///////////////////rotate/////////////////////
		std::vector<Eigen::Vector3d> temp_V;
		temp_V.resize(V.rows());	//temp_V存储旋转后的模型顶点
		for (int i = 0; i < V.rows(); i++) {
			temp_V[i].x() = V.row(i).x();
			temp_V[i].y() = V.row(i).y();
			temp_V[i].z() = V.row(i).z();
		}
		vector<std::vector<Eigen::Vector3d>> temp_new_V;	//temp_new_V存储旋转后的voronoi面边界顶点
		temp_new_V.resize(V.rows());
		for (int i = 0; i < V.rows(); i++) {
			temp_new_V[i].resize(all_voronoi_cells[i].all_points_in_polygon.size());
			for (int k = 0; k < temp_new_V[i].size(); k++) {
				temp_new_V[i][k].x() = all_voronoi_cells[i].all_points_in_polygon[k].x();
				temp_new_V[i][k].y() = all_voronoi_cells[i].all_points_in_polygon[k].y();
				temp_new_V[i][k].z() = all_voronoi_cells[i].all_points_in_polygon[k].z();
			}
		}

		Eigen::Vector3d vectorBefore(0, 0, 1);
		Eigen::Vector3d vectorAfter(sampling_subtractive.sample_points[ori]);
		rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vectorAfter).toRotationMatrix();	//计算从z轴到采样方向的旋转矩阵
		Eigen::Matrix3d rotMatrix_inverse = rotMatrix.inverse();
		for (int i = 0; i < V.rows(); i++)
			temp_V[i] = rotMatrix_inverse * temp_V[i];
		for (int i = 0; i < V.rows(); i++)
			for (int j = 0; j < temp_new_V[i].size(); j++)
				temp_new_V[i][j] = rotMatrix_inverse * temp_new_V[i][j];
		//////////////////////////////////////////////


		/////////////////calculate normal//////////////////
		//vector<vector<Vector3>> all_normal_of_triangles_in_cells;
		vector<Eigen::Vector3d> all_normal_of_cells;	//存储每个voronoi单元的法向量
		//all_normal_of_triangles_in_cells.resize(all_voronoi_cells.size());
		all_normal_of_cells.resize(all_voronoi_cells.size());

		vector<vector<Eigen::Vector3d>> temp_vis;
		for (int i = 0; i < all_voronoi_cells.size(); i++) {
			const auto& cell = all_voronoi_cells[i];
			if (cell.is_available == true) {
				const auto polygon_size = cell.all_points_in_polygon.size();
				const int site_index = cell.site;

				Eigen::Vector3d vn_sum(0.0, 0.0, 0.0);
				for (int j = 0; j < polygon_size; j++) {
					const Eigen::Vector3d v1 = ToVector3(temp_V[site_index]);
					const Eigen::Vector3d v2 = ToVector3(temp_new_V[i][j]);
					const Eigen::Vector3d v3 = ToVector3(temp_new_V[i][(j + 1) % polygon_size]);
					vn_sum += (v2 - v1).cross(v3 - v1);
				}

				all_normal_of_cells[i] = vn_sum;
				if (polygon_size > 0) {
					all_normal_of_cells[i] /= polygon_size;
				}
				all_normal_of_cells[i].normalize();	//voronoi单元法向量归一化

				all_voronoi_cells[i].all_normal_in_all_ori.push_back(all_normal_of_cells[i]);	//将该采样方向下的法向量存入all_voronoi_cells的all_normal_in_all_ori中

				temp_vis.push_back({});
				temp_vis.back().reserve(2);
				Eigen::Vector3d v_site(
					temp_V[site_index].x(),
					temp_V[site_index].y(),
					temp_V[site_index].z());
				Eigen::Vector3d v_normal = v_site + all_normal_of_cells[i] * 3.0;
				temp_vis.back().push_back(v_site);
				temp_vis.back().push_back(v_normal);
			}
		}
		if (ori == 0) {
			Visual vis_normal;//存储法向量的可视化
			//vis_normal.generateModelForRendering_9(temp_vis, file_name);

		}
		//////////////////////////////////////////////////


		///////////////////collision detection////////////////////////
		int cont_accessible_points = 0, cont_unaccessible_points = 0;

		vector<double> max_z_of_cells(all_voronoi_cells.size());	//max_z_of_cells存储每个voronoi单元边界顶点(all_points_in_polygon)的最大z值
		for (int i = 0; i < all_voronoi_cells.size(); i++) {
			max_z_of_cells[i] = MIN_D;
			for (int j = 0; j < all_voronoi_cells[i].all_points_in_polygon.size(); j++)
				max_z_of_cells[i] = max(max_z_of_cells[i], temp_new_V[i][j](2, 0));
		}

		PrepareToolForCollision(the_nozzle);

		for (int i = 0; i < all_voronoi_cells.size(); i++) {	//枚举所有voronoi单元，判断该单元在当前采样方向ori下是否可达
			if (flag_accessible_points[i] == true) {	//如果该单元在之前某个采样方向下已经可达，则跳过
				cont_accessible_points++;
				continue;
			}
			cont_unaccessible_points++;

			if (!all_voronoi_cells[i].is_available) {
				continue;
			}

			// 计算刀尖球中心点坐标
			Eigen::Vector3d center_point = ComputeToolCenter(
				ToVector3(temp_V[all_voronoi_cells[i].site]),
				all_normal_of_cells[i],
				the_nozzle.cylinder_r);

			bool jud_collision = false;
			for (int ii = 0; ii < all_voronoi_cells.size(); ii++) {
				if (i == ii || !all_voronoi_cells[ii].is_available) {
					continue;
				}

				if (CheckToolCollisionWithCell(center_point, temp_new_V[ii], max_z_of_cells[ii], the_nozzle)) {
					jud_collision = true;
					break;
				}
			}

			if (!jud_collision) {
				flag_accessible_points[i] = true;
			}
		}
		num_inaccessible_points = cont_unaccessible_points;
		std::cout << "id of orientation:" << ori << endl;
		std::cout << "number of accessible points:" << cont_accessible_points << endl;
		std::cout << "number of unaccessible points:" << cont_unaccessible_points << endl << endl;
		//break;
	}
	end_time = clock();


	//*******************************************//
	/////////////////find area S//////////////////
	//*******************************************//
	start_time_2 = clock();

	int cont_number_2 = 0;
	std::vector<Eigen::MatrixXd> vis_red_points;	//vis_red_points存储不可达voronoi单元的边界顶点坐标，用于可视化
	std::vector<Eigen::MatrixXd> vis_green_points;	//vis_green_points存储可达voronoi单元的边界顶点坐标，用于可视化
	std::vector<Eigen::Vector3d> temp_V_vis;	//temp_V_vis存储未旋转的模型顶点坐标，用于可视化
	for (int ori = 0; ori < sampling_subtractive.sample_points.size(); ori++) {	//ori枚举所有采样方向
		std::vector<Eigen::Vector3d> temp_V;	//temp_V存储旋转后的模型顶点坐标
		temp_V.resize(V.rows());
		temp_V_vis.resize(V.rows());	//temp_V_vis存储未旋转的模型顶点坐标
		for (int i = 0; i < V.rows(); i++) {
			temp_V[i](0, 0) = V.row(i).x();
			temp_V[i](1, 0) = V.row(i).y();
			temp_V[i](2, 0) = V.row(i).z();
		}
		vector<std::vector<Eigen::Vector3d>> temp_new_V;	//temp_new_V存储旋转后的voronoi面边界顶点坐标
		temp_new_V.resize(V.rows());
		for (int i = 0; i < V.rows(); i++) {
			temp_new_V[i].resize(all_voronoi_cells[i].all_points_in_polygon.size());
			for (int k = 0; k < temp_new_V[i].size(); k++) {
				temp_new_V[i][k](0, 0) = all_voronoi_cells[i].all_points_in_polygon[k].x();
				temp_new_V[i][k](1, 0) = all_voronoi_cells[i].all_points_in_polygon[k].y();
				temp_new_V[i][k](2, 0) = all_voronoi_cells[i].all_points_in_polygon[k].z();
			}
		}
		temp_V_vis = temp_V;

		Eigen::Vector3d vectorBefore(0, 0, 1);
		Eigen::Vector3d vectorAfter(sampling_subtractive.sample_points[ori]);
		rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vectorAfter).toRotationMatrix();
		for (int i = 0; i < V.rows(); i++)
			temp_V[i] = rotMatrix.inverse() * temp_V[i];	//将模型顶点旋转到以采样方向为z轴的坐标系下
		for (int i = 0; i < V.rows(); i++)
			for (int j = 0; j < temp_new_V[i].size(); j++)
				temp_new_V[i][j] = rotMatrix.inverse() * temp_new_V[i][j];	//将voronoi面边界顶点旋转到以采样方向为z轴的坐标系下

		//cout << "b" << endl;
		/////////////////calculate normal//////////////////
		//vector<vector<Vector3>> all_normal_of_triangles_in_cells;
		vector<Eigen::Vector3d> all_normal_of_cells;	//仅计算不可达voronoi单元的法向量
		//all_normal_of_triangles_in_cells.resize(all_voronoi_cells.size());
		all_normal_of_cells.resize(all_voronoi_cells.size());
		for (int i = 0; i < all_voronoi_cells.size(); i++) {
			if (all_voronoi_cells[i].is_available == true && flag_accessible_points[i] == false) {	//仅计算不可达voronoi单元的法向量
				//all_normal_of_triangles_in_cells[i].resize(all_voronoi_cells[i].all_points_in_polygon.size());
				//all_normal_of_cells[i].x() = all_normal_of_cells[i].y() = all_normal_of_cells[i].z() = 0;
				//for (int j = 0; j < 1; j++) {
				//	Eigen::Vector3d v1 = ToVector3(temp_V[all_voronoi_cells[i].site]);
				//	Eigen::Vector3d v2 = ToVector3(temp_new_V[i][j]);
				//	Eigen::Vector3d v3 = ToVector3(temp_new_V[i][(j + 1) % all_voronoi_cells[i].all_points_in_polygon.size()]);
				//	Eigen::Vector3d vn = (v2 - v1).cross(v3 - v1);
				//	//all_normal_of_triangles_in_cells[i][j] = vn;
				//	all_normal_of_cells[i].x() += vn.x();
				//	all_normal_of_cells[i].y() += vn.y();
				//	all_normal_of_cells[i].z() += vn.z();
				//}
				//all_normal_of_cells[i] /= all_voronoi_cells[i].all_points_in_polygon.size();
				//all_normal_of_cells[i].normalize();
				all_normal_of_cells[i] = all_voronoi_cells[i].all_normal_in_all_ori[ori];	//直接使用之前计算好的法向量
			}
		}
		//////////////////////////////////////////////////
		//cout << "c" << endl;
		//ofstream vis_unaccessible_points(file_name + "_unaccessible_points.obj");

		file_name_2 = file_name + to_string(cont_num);
		vector<double> max_z_of_cells(all_voronoi_cells.size());	//max_z_of_cells存储每个voronoi单元边界顶点(all_points_in_polygon)的最大z值
		for (int i = 0; i < all_voronoi_cells.size(); i++) {
			max_z_of_cells[i] = MIN_D;
			for (int j = 0; j < all_voronoi_cells[i].all_points_in_polygon.size(); j++)
				max_z_of_cells[i] = max(max_z_of_cells[i], temp_new_V[i][j](2, 0));
		}

		int cont_unaccessible_points = 0;

		PrepareToolForCollision(the_nozzle);

		for (int i = 0; i < all_voronoi_cells.size(); i++) {
			if (flag_accessible_points[i] == true) {	//跳过可达voronoi单元
				continue;
			}

			if (!all_voronoi_cells[i].is_available) {
				continue;
			}

			// 计算不可达voronoi单元i的刀尖球中心点坐标
			Eigen::Vector3d center_point = ComputeToolCenter(
				ToVector3(temp_V[all_voronoi_cells[i].site]),
				all_normal_of_cells[i],
				the_nozzle.cylinder_r);

			int index_insert = 0;
			if (ori == 0) {	//在首个方向 ori==0 时，初始化该不可达单元的 area_S 容器与索引映射
				vector<area_S> temp_vec_area_s;
				//temp_vec_area_s.resize(1);
				all_the_area_S.push_back(temp_vec_area_s);
				map_S_and_vertex.insert({ cont_unaccessible_points,i });
				map_S_and_vertex_inv.insert({ i, cont_unaccessible_points });
				cont_unaccessible_points++;	//记录不可达单元数量,第(cont_unaccessible_points-1)个不可达voronoi单元的索引为i
			}
			else {
				index_insert = map_S_and_vertex_inv[i];
			}

			//vis_unaccessible_points << "v " << temp_V[i](0,0) << " " << temp_V[i](1, 0) << " " << temp_V[i](2, 0) << endl;
			start_time_4 = clock();
			if (ori == 0)
				vis_red_points.push_back(temp_V_vis[i]);	//将不可达voronoi单元i的网格顶点存入vis_red_points，用于可视化
			end_time_4 = clock();
			sum_time += (double(end_time_4) - double(start_time_4)) / CLOCKS_PER_SEC;

			for (int ii = 0; ii < all_voronoi_cells.size(); ii++) {	//枚举所有voronoi单元，判断刀尖放在当前单元i、刀具方向是ori方向时，是否与其他单元发生碰撞
				if (i == ii || !all_voronoi_cells[ii].is_available) {
					continue;
				}

				if (CheckToolCollisionWithCell(center_point, temp_new_V[ii], max_z_of_cells[ii], the_nozzle)) {
					area_S temp_area_S(ii, ori);
					if (ori == 0) {
						all_the_area_S.back().push_back(temp_area_S);
					}
					else {
						all_the_area_S[index_insert].push_back(temp_area_S);
						cont_number_2++;
					}
					test_cont_num++;
				}
			}
		}

		//cout << "d" << endl;
		Visual Vis;
		//Vis.generateModelForRendering(temp_layers, file_name_2);
		//Vis.generateModelForRendering_3(vectorAfter, file_name_2, vis_points);
		//break;

		cont_num++;
	}
	//*******************************************//
	//cout << endl << all_the_area_S.size() << endl;
	//cout << endl << cont_number_2 << endl;
	end_time_2 = clock();

	////////////get covering points(green points)////////////
	//only consider covering for area S,save the index of all_the_area_S
	start_time_3 = clock();
	cout << test_cont_num << "&&&&" << endl;
	int cont_covering_points = 0;
	for (int i = 0; i < all_the_area_S.size(); i++) {	//枚举第i个不可达voronoi单元的area_S
		for (int j = 0; j < all_the_area_S[i].size(); j++) {	//枚举里面的所有area_S元素
			if (flag_covering_points[all_the_area_S[i][j].pointId] == false) {	//如果该voronoi单元还没有在flag_covering_points中被标记
				map_covering_points_and_vertex.insert({ cont_covering_points,all_the_area_S[i][j].pointId });	//建立被碰撞点与voronoi单元索引的映射
				map_covering_points_and_vertex_inv.insert({ all_the_area_S[i][j].pointId, cont_covering_points });
				vector<area_S> temp_vec_area_s;
				all_the_covering_points.push_back(temp_vec_area_s);	//在all_the_covering_points中添加一个新的容器

				area_S temp_covering_point(i, all_the_area_S[i][j].oriId);	//
				all_the_covering_points[all_the_covering_points.size() - 1].push_back(temp_covering_point);
				flag_covering_points[all_the_area_S[i][j].pointId] = true;
				cont_covering_points++;
			}
			else {
				int index_insert = map_covering_points_and_vertex_inv[all_the_area_S[i][j].pointId];
				area_S temp_covering_point(i, all_the_area_S[i][j].oriId);
				all_the_covering_points[index_insert].push_back(temp_covering_point);
			}
		}
	}
	int max_size = -MAX_I;
	for (int i = 0; i < all_the_covering_points.size(); i++) {
		max_size = max(max_size, int(all_the_covering_points[i].size()));
		vis_green_points.push_back(temp_V_vis[map_covering_points_and_vertex[i]]);	//将覆盖点对应的voronoi单元顶点存入vis_green_points，用于可视化
	}
	vector<double> color_map;
	color_map.resize(all_the_covering_points.size());
	for (int i = 0; i < all_the_covering_points.size(); i++)
		color_map[i] = double(all_the_covering_points[i].size()) / double(max_size);	//归一化颜色值，便于可视化，被更多个area S覆盖的点颜色更亮
	end_time_3 = clock();
	/*vector<vector<area_S>> temp_all_the_covering_points = all_the_covering_points;
	map<int, int> temp_map_covering_points_and_vertex = map_covering_points_and_vertex;
	for(int i = 0;i< temp_all_the_covering_points.size();i++)
		for (int j = i + 1; j < temp_all_the_covering_points.size(); j++)
			if (temp_all_the_covering_points[i].size() < temp_all_the_covering_points[j].size()) {
				swap(temp_all_the_covering_points[i], temp_all_the_covering_points[j]);
				swap(temp_map_covering_points_and_vertex[i], temp_map_covering_points_and_vertex[j]);
			}

	for(int i =0;i < temp_all_the_covering_points.size() / 10;i++)
		vis_green_points.push_back(temp_V_vis[temp_map_covering_points_and_vertex[i]]);*/
		/////////////////////////////////////////////////////////
	const std::string accessibility_point_vis_prefix =
		UnaccessibleVisPath(MakeAccessibilityDebugFileToken(file_name));
	if (open_vis_red_points == true)
		createRedBalls(accessibility_point_vis_prefix, vis_red_points);
	if (open_vis_green_points == true)
		createGreenBalls(accessibility_point_vis_prefix, vis_green_points, color_map);
	accessibility_visualization::WriteInitialSubtractiveAllDirectionDebugVisualization(
		UnaccessibleVisDir(),
		MakeAccessibilityDebugFileToken(file_name),
		std::max(0.25, dh * 0.15),
		V,
		map_S_and_vertex);
	//accessibility_visualization::WriteSubtractiveAllDirectionReasonDiagnostics(
	//	UnaccessibleVisDir(),
	//	MakeAccessibilityDebugFileToken(file_name),
	//	V,
	//	all_voronoi_cells,
	//	sampling_subtractive.sample_points,
	//	map_S_and_vertex,
	//	the_nozzle);

	std::cout << "Collision Detection For Subtractive Manufacturing have done" << endl;
	/*std::cout << "&&&time&&& Collision detection: " << double(end_time - start_time) / CLOCKS_PER_SEC << std::endl;
	std::cout << "&&&time&&& Subtractive dependency graph: " << double(end_time_2 - start_time_2) / CLOCKS_PER_SEC << std::endl;
	std::cout << "&&&time&&& Green points generation: " << double(end_time_3 - start_time_3) / CLOCKS_PER_SEC << std::endl;*/
	time_build_subtractive_graph = double(end_time - start_time) / CLOCKS_PER_SEC + double(end_time_2 - start_time_2) / CLOCKS_PER_SEC + double(end_time_3 - start_time_3) / CLOCKS_PER_SEC;
	std::cout << "&&&time&&& Build subtractive graph: " << double(end_time - start_time) / CLOCKS_PER_SEC + double(end_time_2 - start_time_2) / CLOCKS_PER_SEC + double(end_time_3 - start_time_3) / CLOCKS_PER_SEC << std::endl;

	if (all_the_area_S.size() == 0) {
		cout << file_name << " 无不可达点！！！" << endl;

		Eigen::MatrixXd V;
		Eigen::MatrixXi F;
		Eigen::MatrixXi N;

		igl::readOBJ(file_name + ".obj", V, F);
		igl::writeOBJ(file_name + "-1_0_current.obj", V, F);
		return -1;
	}
	return 0;
}

void HybridManufacturing::GetALLFragileVertex(SAMPLE_ON_BALL sampling)
{
	for (int ori = 0; ori < sampling.sample_points.size(); ori++) {
		//rotating the blocks and then slicing//
		vector<Eigen::Vector3d> temp_V;
		temp_V.resize(V.rows());
		for (int i = 0; i < V.rows(); i++) {
			temp_V[i](0, 0) = V.row(i).x();
			temp_V[i](1, 0) = V.row(i).y();
			temp_V[i](2, 0) = V.row(i).z();
		}
		Eigen::Vector3d vectorBefore(0, 0, 1);
		Eigen::Vector3d vectorAfter(sampling.sample_points[ori]);
		Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vectorAfter).toRotationMatrix();
		for (int i = 0; i < V.rows(); i++)
			temp_V[i] = rotMatrix.inverse() * temp_V[i];

		vector<Eigen::Vector3d> normal_V;
		vector<bool> temp_fragile_V;
		normal_V.resize(V.rows());
		temp_fragile_V.resize(V.rows());
		for (int i = 0; i < V.rows(); i++) {
			normal_V[i].x() = normal_V[i].y() = normal_V[i].z() = 0;
		}
		for (int i = 0; i < F.rows(); i++) {
			VEctor v1, v2, v3;
			v1[0] = temp_V[F(i, 0)](0, 0);
			v1[1] = temp_V[F(i, 0)](1, 0);
			v1[2] = temp_V[F(i, 0)](2, 0);
			v2[0] = temp_V[F(i, 1)](0, 0);
			v2[1] = temp_V[F(i, 1)](1, 0);
			v2[2] = temp_V[F(i, 1)](2, 0);
			v3[0] = temp_V[F(i, 2)](0, 0);
			v3[1] = temp_V[F(i, 2)](1, 0);
			v3[2] = temp_V[F(i, 2)](2, 0);
			Eigen::Vector3d temp_normal;
			temp_normal = faceNormal(v1, v2, v3);
			temp_normal.normalize();
			normal_V[F(i, 0)] += temp_normal;
			normal_V[F(i, 1)] += temp_normal;
			normal_V[F(i, 2)] += temp_normal;
		}
		for (int i = 0; i < V.rows(); i++) {
			normal_V[i].normalize();
			if (normal_V[i].z() > 0.95)
				temp_fragile_V[i] = true;
			else
				temp_fragile_V[i] = false;
		}
		is_fragile_V.push_back(temp_fragile_V);
	}
}

void HybridManufacturing::detect_collision_with_printing_platform(
	int& index,
	vector<int>& candidate_nodes,
	OrientationScores& all_calculated_value,
	CutLayerVector all_cut_layers,
	Eigen::Vector3d ori_now,
	nozzle the_nozzle)
{

	Eigen::Vector3d center_point(0.0, 0.0, 0.0);
	for (int i = 0; i < V_bottom.size(); i++) {
		center_point += V_bottom[i];
	}
	center_point /= V_bottom.size();
	double circle_r = 160;

	if (ori_now.x() == 0 && ori_now.y() == 0 && ori_now.z() == 1)
		return;
	//if (all_calculated_value[index].number_of_remaining_face < terminate_threshold_of_number_of_faces) 
		//return;
	Eigen::Vector3d vectorBefore(0, 0, 1);
	Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(ori_now, vectorBefore).toRotationMatrix();
	vector<vector<Eigen::Vector3d>> rotate_all_cut_layers;
	rotate_all_cut_layers.resize(all_cut_layers.size());
	for (int i = 0; i < all_cut_layers.size(); i++) {
		rotate_all_cut_layers[i].resize(all_cut_layers[i].size());
		for (int j = 0; j < all_cut_layers[i].size(); j++) {
			rotate_all_cut_layers[i][j](0, 0) = all_cut_layers[i][j].x;
			rotate_all_cut_layers[i][j](1, 0) = all_cut_layers[i][j].y;
			rotate_all_cut_layers[i][j](2, 0) = all_cut_layers[i][j].z;
			rotate_all_cut_layers[i][j] = rotMatrix.inverse() * rotate_all_cut_layers[i][j];
		}
	}

	for (int i = 0; i < rotate_all_cut_layers.size(); i++) {
		for (int j = 0; j < rotate_all_cut_layers[i].size(); j++) {
			if ((abs(rotate_all_cut_layers[i][j](2, 0) - center_point.z()) < the_nozzle.lowwer_surface_r) && (pow(rotate_all_cut_layers[i][j](0, 0) - center_point.x(), 2) + pow(rotate_all_cut_layers[i][j](1, 0) - center_point.y(), 2) - pow(circle_r, 2) < 0)) {
				candidate_nodes.erase(candidate_nodes.begin() + index);
				all_calculated_value.erase(all_calculated_value.begin() + index);
				index--;
				return;
			}
		}
	}
}

all_value HybridManufacturing::GainMesh(
	Slicer_2& slicer,
	CutLayerVector all_cut_layers,
	Eigen::Vector3d vector_after,
	int height_of_beam_search,
	int cont_number_of_queue,
	int index_of_pre_node,
	vector<int> all_cut_layers_dependency_layer,
	bool flag_is_continue_block,
	int id_continue)
{
	clock_t start_time, end_time;

	double sum_area = 0;
	double sum_value_of_manufacturing_require = 0;
	double sum_not_self_support_area = 0;
	//double area_threshold = 5.0;
	double ratio_of_area_threshold = 0.05;  //0.015 //0.15  //0.08   //0.03  //0.05
	double max_self_support_slope_angle;
	//cout << "#####" << height_of_beam_search << " " << id_continue << endl;
	//if(height_of_beam_search != 7 || id_continue != 1)s
	//	max_self_support_slope_angle = PI / 10;  //PI/3.6  
	//else
	max_self_support_slope_angle = PI / 3.6;  //PI/3.6  
	double sum_area_threshold = 1000;  //200 //1000
	if (vector_after.x() == 0 && vector_after.y() == 0) {
		sum_area_threshold = 100000;
		ratio_of_area_threshold = 0.1;
	}

	/*if(height_of_beam_search == 4)
		ratio_of_area_threshold = 0.2;*/

	std::vector<bool> jud_triangle_have_been_added;
	vector<vector<int>> cutting_plane_points;
	cutting_plane_points.resize(all_cut_layers.size());
	//////sort cut layers//////
	for (int i = 0; i < all_cut_layers.size(); i++)
		for (int j = i + 1; j < all_cut_layers.size(); j++) {

			if (all_cut_layers[i][0].z > all_cut_layers[j][0].z) {
				swap(all_cut_layers[i], all_cut_layers[j]);
				swap(all_cut_layers_dependency_layer[i], all_cut_layers_dependency_layer[j]);
			}
		}
	///////////////////////////


	//need rotate first//
	Eigen::Vector3d vectorBefore(0, 0, 1);
	Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vector_after).toRotationMatrix();
	for (int i = 0; i < slicer.positions.size(); i++) {
		Eigen::Vector3d temp_V;
		temp_V(0, 0) = slicer.positions[i][0];
		temp_V(1, 0) = slicer.positions[i][1];
		temp_V(2, 0) = slicer.positions[i][2];
		temp_V = rotMatrix.inverse() * temp_V;
		slicer.positions[i][0] = temp_V(0, 0);
		slicer.positions[i][1] = temp_V(1, 0);
		slicer.positions[i][2] = temp_V(2, 0);
	}
	/////////////////////


	//--------------------cut-----------------------//
	slicer.normal[0] = 0;
	slicer.normal[1] = 0;
	slicer.normal[2] = 1;
	//cout << endl << "pp" << endl;
	for (int i = 0; i < all_cut_layers.size(); i++) {
		slicer.origin[0] = 0;
		slicer.origin[1] = 0;
		slicer.origin[2] = all_cut_layers[i][0].z;
		slicer.cut();
	}
	//cout << endl << "qq" << endl;
	std::vector<TRiangle> ori_triangle = slicer.triangles;
	jud_triangle_have_been_added.resize(slicer.triangles.size());
	for (int i = 0; i < slicer.triangles.size(); i++)
		jud_triangle_have_been_added[i] = false;

	Slicer_2 all_slicer;

	all_slicer.positions = slicer.positions;
	all_slicer.triangles = ori_triangle;


	//--------------------save candidate_triangles-----------------------//
	int current_index = 0;
	std::vector<VEctor> min_z_point;
	std::vector<double> min_z_triangle;
	vector<int> index_of_min_point_in_triangle;
	std::vector<TRiangle> candidate_triangles;
	std::vector<int> id_triangles;
	std::vector<TRiangle> remove_triangles;
	candidate_triangles.clear();
	id_triangles.clear();
	double min_z_all_cut_layers = 9999999;
	for (int t = 0; t < all_cut_layers.size(); t++)
		min_z_all_cut_layers = min(min_z_all_cut_layers, all_cut_layers[t][0].z);
	for (int i = 0; i < all_slicer.triangles.size(); i++) {
		double temp_min_z_triangle = 9999999;
		for (int k = 0; k < 3; k++) {
			temp_min_z_triangle = min(all_slicer.positions[all_slicer.triangles[i][k]][2], temp_min_z_triangle);
		}
		if (temp_min_z_triangle + 0.001 >= min_z_all_cut_layers) {
			candidate_triangles.push_back(all_slicer.triangles[i]);
			id_triangles.push_back(i);
		}
	}
	for (int i = 0; i < candidate_triangles.size(); i++) {
		double temp_min_z_triangle = 9999999;
		VEctor temp_min_z_point = all_slicer.positions[candidate_triangles[i][0]];
		int temp_index_of_min_point_in_triangle;
		for (int k = 0; k < 3; k++) {
			if (all_slicer.positions[candidate_triangles[i][k]][2] < temp_min_z_triangle) {
				temp_min_z_triangle = all_slicer.positions[candidate_triangles[i][k]][2];
				temp_min_z_point = all_slicer.positions[candidate_triangles[i][k]];
				temp_index_of_min_point_in_triangle = candidate_triangles[i][k];
			}
		}
		min_z_triangle.push_back(temp_min_z_triangle);
		min_z_point.push_back(temp_min_z_point);
		index_of_min_point_in_triangle.push_back(temp_index_of_min_point_in_triangle);
	}

	quick_sort(candidate_triangles, 0, candidate_triangles.size() - 1, id_triangles, min_z_triangle, min_z_point);


	////--------------------save OPP_triangles one by one-----------------------//
	vector<vector<TRiangle>> all_furcation_of_blocks;
	all_furcation_of_blocks.resize(all_cut_layers.size());
	vector<int> id_of_furcation_of_blocks;
	vector<int> save_current_index;
	for (int t = 0; t < all_cut_layers.size(); t++) {
		current_index = 0;   //该方式也许比较慢
		double boundary_bottom = 999999, boundary_left = 999999, boundary_top = -999999, boundary_right = -999999;
		for (int i = 0; i < all_cut_layers[t].size(); i++) {
			boundary_top = std::max(boundary_top, all_cut_layers[t][i].y);
			boundary_bottom = std::min(boundary_bottom, all_cut_layers[t][i].y);
			boundary_right = std::max(boundary_right, all_cut_layers[t][i].x);
			boundary_left = std::min(boundary_left, all_cut_layers[t][i].x);
		}
		Eigen::Vector2d current_triangle_point;
		Eigen::Vector2d current_layer_point;
		for (; current_index < candidate_triangles.size(); current_index++) {
			if (abs(min_z_triangle[current_index] - all_cut_layers[t][0].z) > 0.001) {
				if (min_z_triangle[current_index] > all_cut_layers[t][0].z)
					break;
			}
			else {
				bool jud_is_boundary_point = false;
				//for (int i = 0; i < all_cut_layers[t].size(); i++) {
				int cont_inside_boundary = 0;
				for (int k = 0; k < 3; k++) {
					if (all_slicer.positions[candidate_triangles[current_index][k]][0] + 0.01 >= boundary_left && all_slicer.positions[candidate_triangles[current_index][k]][0] - 0.01 <= boundary_right
						&& all_slicer.positions[candidate_triangles[current_index][k]][1] + 0.01 >= boundary_bottom && all_slicer.positions[candidate_triangles[current_index][k]][1] - 0.01 <= boundary_top) {
						cont_inside_boundary++;
					}
				}
				if (cont_inside_boundary >= 2)
					jud_is_boundary_point = true;
				/*if (min_z_point[current_index][0] + 0.01 >= boundary_left && min_z_point[current_index][0] - 0.01 <= boundary_right &&
					min_z_point[current_index][1] + 0.01 >= boundary_bottom && min_z_point[current_index][1] - 0.01 <= boundary_top) {
					jud_is_boundary_point = true;
					break;
				}*/
				//}
				if (jud_is_boundary_point == true) {
					for (int j = 0; j < all_cut_layers[t].size(); j++) {
						current_layer_point.x() = all_cut_layers[t][j].x;
						current_layer_point.y() = all_cut_layers[t][j].y;
						current_triangle_point.x() = min_z_point[current_index][0];
						current_triangle_point.y() = min_z_point[current_index][1];
						if ((current_layer_point - current_triangle_point).norm() < 0.2) {  //4.0
							//if (height_of_beam_search != 2)
							remove_triangles.push_back(candidate_triangles[current_index]);
							id_of_furcation_of_blocks.push_back(t);
							all_furcation_of_blocks[t].push_back(candidate_triangles[current_index]);

							Eigen::Vector3d face_normal = faceNormal(all_slicer.positions[remove_triangles[remove_triangles.size() - 1][0]], all_slicer.positions[remove_triangles[remove_triangles.size() - 1][1]], all_slicer.positions[remove_triangles[remove_triangles.size() - 1][2]]);
							face_normal.normalize();
							Eigen::Vector3d base_normal(0, 0, 1);
							int jud_self_support = (face_normal.dot(base_normal) + sin(max_self_support_slope_angle) >= 0);
							double current_area = triangleArea(all_slicer.positions[remove_triangles[remove_triangles.size() - 1][0]], all_slicer.positions[remove_triangles[remove_triangles.size() - 1][1]], all_slicer.positions[remove_triangles[remove_triangles.size() - 1][2]]);
							sum_area += current_area;
							//sum_value_of_manufacturing_require += current_area * calculate_manufacturing_require(all_slicer, remove_triangles);
							if (jud_self_support == false) {
								sum_not_self_support_area += current_area;
								if (sum_not_self_support_area / sum_area > ratio_of_area_threshold && sum_area >= sum_area_threshold) {
									all_value temp_all_value;
									temp_all_value.value_of_self_support = 0;
									return temp_all_value;
								}
							}

							save_current_index.push_back(current_index);
							jud_triangle_have_been_added[id_triangles[current_index]] = true;
							break;
						}
					}
				}
			}

		}
	}
	current_index = 0;  //current_index = 0;
	//std::cout << "**" << current_index << " " << candidate_triangles.size() << " " << remove_triangles.size() << endl;

	start_time = clock();
	while (1) {
		bool flag_break = false;
		for (int i = 0; i < remove_triangles.size(); i++)
		{
			if (jud_triangle_have_been_added[id_triangles[current_index]] == false && min_z_triangle[current_index] >= min_z_triangle[save_current_index[i]]) {
				for (int j = 0; j < 3; j++) {
					for (int k = 0; k < 3; k++) {
						if (candidate_triangles[current_index][j] == remove_triangles[i][k]) {
							//if(height_of_beam_search != 2)
							remove_triangles.push_back(candidate_triangles[current_index]);
							id_of_furcation_of_blocks.push_back(id_of_furcation_of_blocks[i]);
							all_furcation_of_blocks[id_of_furcation_of_blocks[i]].push_back(candidate_triangles[current_index]);
							//////all_furcation_of_blocks[t].push_back(candidate_triangles[current_index]);

							Eigen::Vector3d face_normal = faceNormal(all_slicer.positions[remove_triangles[remove_triangles.size() - 1][0]], all_slicer.positions[remove_triangles[remove_triangles.size() - 1][1]], all_slicer.positions[remove_triangles[remove_triangles.size() - 1][2]]);
							face_normal.normalize();
							Eigen::Vector3d base_normal(0, 0, 1);
							int jud_self_support = (face_normal.dot(base_normal) + sin(max_self_support_slope_angle) >= 0);
							double current_area = triangleArea(all_slicer.positions[remove_triangles[remove_triangles.size() - 1][0]], all_slicer.positions[remove_triangles[remove_triangles.size() - 1][1]], all_slicer.positions[remove_triangles[remove_triangles.size() - 1][2]]);
							sum_area += current_area;

							if (jud_self_support == false) {
								sum_not_self_support_area += current_area;
								if (sum_not_self_support_area / sum_area > ratio_of_area_threshold && sum_area >= sum_area_threshold) {
									all_value temp_all_value;
									temp_all_value.value_of_self_support = 0;
									return temp_all_value;
								}
							}

							save_current_index.push_back(current_index);
							jud_triangle_have_been_added[id_triangles[current_index]] = true;
							flag_break = true;
							break;
						}
					}
					if (flag_break == true)
						break;
				}
			}
			if (flag_break == true)
				break;
		}
		current_index++;
		if (current_index >= candidate_triangles.size())
			break;
	}end_time = clock();
	//add remaining face, set a distance threshold Dis, only a face less Dis from other cut layers and no other dependency layer exist, the face is consider to remaining face
	/*double Dis = dh * 4;
	for (int i = 0; i < jud_triangle_have_been_added.size(); i++) {
		if (jud_triangle_have_been_added[i] == false) {
			int id_layer;
			double min_dis = 99999999;
			for (int t = 0; t < all_cut_layers.size(); t++) {
				if (all_cut_layers[t][0].z - slicer.positions[slicer.triangles[i][0]][2] <= 2 * dh && all_cut_layers[t][0].z - slicer.positions[slicer.triangles[i][0]][2] >= -0.001) {
					for (int j = 0; j < all_cut_layers[t].size(); j++) {
						cv::Point3d current_triangle_point(slicer.positions[slicer.triangles[i][0]][0], slicer.positions[slicer.triangles[i][0]][1], slicer.positions[slicer.triangles[i][0]][2]);
						cv::Point3d current_layer_point(all_cut_layers[t][j].x, all_cut_layers[t][j].y, all_cut_layers[t][0].z);

						double distance = distance_3d(current_triangle_point, current_layer_point);
						if (distance < min_dis) {
							min_dis = distance;
							id_layer = t;
						}
					}
				}
			}
			if (min_dis < Dis && all_cut_layers_dependency_layer[id_layer] == 0)
				remove_triangles.push_back(slicer.triangles[i]);
		}
	}*/

	all_value all_calculated_value;
	all_calculated_value.number_of_remaining_face = all_slicer.triangles.size() - remove_triangles.size();
	all_slicer.triangles = remove_triangles;


	all_calculated_value.value_of_projected = calculate_projected_area(all_slicer, all_furcation_of_blocks, all_cut_layers);

	double min_z = 9999999, max_z = -9999999;
	double value_of_manufacturing_require;
	for (int i = 0; i < all_slicer.positions.size(); i++) {
		min_z = min(min_z, all_slicer.positions[i][2]);
		max_z = max(max_z, all_slicer.positions[i][2]);
	}
	for (int i = 0; i < all_slicer.triangles.size(); i++) {
		double current_area = triangleArea(all_slicer.positions[all_slicer.triangles[i][0]], all_slicer.positions[all_slicer.triangles[i][1]], all_slicer.positions[all_slicer.triangles[i][2]]);
		value_of_manufacturing_require = vasco::calculateManufacturingRequire(all_slicer, i, min_z, max_z);
		sum_value_of_manufacturing_require += value_of_manufacturing_require * current_area;
	}

	all_calculated_value.value_of_self_support = 1 - sum_not_self_support_area / sum_area;
	all_calculated_value.large_base = sum_value_of_manufacturing_require / sum_area;


	//将slicer旋转回原始位置
	RotateSlicerPositions(all_slicer, vector_after, vectorBefore);

	//auto filtered_patch = FilterSurfaceRemoveTriangles(all_slicer, remove_triangles);
	//if (!all_slicer.positions.empty() && !filtered_patch.empty()) {
	//	const std::string mesh_name = "all_slicer_" + std::to_string(height_of_beam_search) + "_" + std::to_string(cont_number_of_queue);
	//	polyscope::registerSurfaceMesh(mesh_name, all_slicer.positions, filtered_patch);
	//}

	//if (height_of_beam_search > 1) {
	//	subtractive_remove_output(filtered_patch, all_slicer, height_of_beam_search);
	//	cutter cutter1;
	//	int sub_patches = subtractive_accessibility_decomposition_local(height_of_beam_search, cutter1);
	//	std::cout << "[GainMesh] sub_patches: " << sub_patches <<
	//		" cont_number_of_queue " << cont_number_of_queue <<
	//		" index_of_pre_node " << index_of_pre_node <<
	//		" height_of_beam_search " << height_of_beam_search << endl;
	//	all_calculated_value.value_of_sub_patches = 10.0 / sub_patches;
	//	//polyscope::show();

	//	//cout << "()()()(" << double(end_time - start_time) / CLOCKS_PER_SEC << endl;
	//}
	return all_calculated_value;
}

void HybridManufacturing::SortCutLayersByHeight(
	vector<vector<cv::Point3d>>& all_cut_layers,
	vector<int>& all_cut_layers_dependency_layer,
	vector<int>& flag_cut_layers_is_hole,
	map<int, int>& follow_index) const
{
	follow_index.clear();
	for (int i = 0; i < all_cut_layers.size(); i++)
		follow_index.insert({ i,i });	//初始化follow_index
	for (int i = 0; i < all_cut_layers.size(); i++) {
		for (int j = i + 1; j < all_cut_layers.size(); j++) {	//按照z值从小到大排序cut layers
			if (all_cut_layers[i][0].z > all_cut_layers[j][0].z) {
				swap(all_cut_layers[i], all_cut_layers[j]);	//交换cut layers
				swap(follow_index[i], follow_index[j]);	//交换对应关系
				swap(flag_cut_layers_is_hole[i], flag_cut_layers_is_hole[j]);	//交换是否为孔洞的标记
				swap(all_cut_layers_dependency_layer[i], all_cut_layers_dependency_layer[j]);	//交换对应的依赖层数量
			}
		}
	}
}

Slicer_2 HybridManufacturing::LoadSlicerForCutMesh(
	bool flag_is_continue_block,
	int height_of_beam_search,
	int index_of_pre_node,
	int id_continue) const
{
	Slicer_2 slicer;
	if (!flag_is_continue_block) {
		slicer.load(file_name + "-" + to_string(height_of_beam_search - 1) + "_" + to_string(index_of_pre_node) + ".obj");
		cout << "&" << endl;
	}
	else {
		slicer.load(file_name + "-" + to_string(height_of_beam_search - 1) + "_" + to_string(index_of_pre_node) + "_" + to_string(id_continue - 1) + "_subblock.obj");
		cout << "*" << endl;
	}
	return slicer;
}

void HybridManufacturing::RotateSlicerPositions(
	Slicer_2& slicer,
	const Eigen::Vector3d& vector_before,
	const Eigen::Vector3d& vector_after) const
{
	const Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vector_before, vector_after).toRotationMatrix();
	for (int i = 0; i < slicer.positions.size(); i++) {
		Eigen::MatrixXd temp_V(3, 1);
		temp_V(0, 0) = slicer.positions[i][0];
		temp_V(1, 0) = slicer.positions[i][1];
		temp_V(2, 0) = slicer.positions[i][2];
		temp_V = rotMatrix.inverse() * temp_V;
		slicer.positions[i][0] = temp_V(0, 0);
		slicer.positions[i][1] = temp_V(1, 0);
		slicer.positions[i][2] = temp_V(2, 0);
	}
}

bool HybridManufacturing::SaveRemovedSolidAfterCut(
	const Slicer_2& remaining_slicer,
	const std::vector<TRiangle>& remove_triangles,
	const std::vector<int>& contact_face_ids,
	const Eigen::Vector3d& vector_after,
	const Eigen::Vector3d& vector_before,
	int height_of_beam_search,
	int cont_number_of_queue,
	bool judge_continue_additive,
	int id_continue) const
{
	if (remove_triangles.empty()) {
		std::cout << "[HybridManufacturing::CutMesh] Skip removed solid output because no triangles were removed."
			<< std::endl;
		return false;
	}

	Slicer_2 removed_solid;
	removed_solid.positions = remaining_slicer.positions;
	removed_solid.triangles = remove_triangles;
	removed_solid.triangles.reserve(remove_triangles.size() + contact_face_ids.size());

	std::size_t added_contact_face_count = 0;
	std::size_t invalid_contact_face_count = 0;
	for (int face_id : contact_face_ids) {
		if (face_id < 0 || face_id >= static_cast<int>(remaining_slicer.triangles.size())) {
			++invalid_contact_face_count;
			continue;
		}

		TRiangle contact_face = remaining_slicer.triangles[face_id];
		std::swap(contact_face[1], contact_face[2]);
		removed_solid.triangles.push_back(contact_face);
		++added_contact_face_count;
	}

	RotateSlicerPositions(removed_solid, vector_after, vector_before);

	std::string output_file = file_name + "-"
		+ to_string(height_of_beam_search) + "_"
		+ to_string(cont_number_of_queue);
	if (judge_continue_additive) {
		output_file += "_" + to_string(id_continue) + "_removed_solid_subblock.obj";
	}
	else {
		output_file += "_removed_solid.obj";
	}

	if (!removed_solid.save(output_file)) {
		std::cout << "[HybridManufacturing::CutMesh] Failed to save removed solid: "
			<< output_file << std::endl;
		return false;
	}

	std::cout << "[HybridManufacturing::CutMesh] Removed solid saved: "
		<< output_file
		<< ", surface_faces=" << remove_triangles.size()
		<< ", reversed_contact_faces=" << added_contact_face_count
		<< ", invalid_contact_face_ids=" << invalid_contact_face_count
		<< std::endl;
	return true;
}

void HybridManufacturing::AppendCutBoundaryCaps(
	std::vector<std::vector<std::vector<int>>> loops_by_layer,
	const std::vector<int>& flag_cut_layers_is_hole,
	const std::map<int, int>& follow_index,
	const std::vector<int>& all_cut_layers_dependency_layer,
	Slicer_2& all_slicer,
	std::vector<int>& all_cap_face_ids,
	std::vector<int>& id_contact_faces) const
{
	std::map<vasco::contact_triangulation::QuantizedPointKey, int> point_index_map;
	for (const auto& tri : all_slicer.triangles) {
		for (int k = 0; k < 3; ++k) {
			point_index_map.emplace(
				vasco::contact_triangulation::MakeQuantizedKey(Point_3(
					all_slicer.positions[tri[k]][0],
					all_slicer.positions[tri[k]][1],
					all_slicer.positions[tri[k]][2])),
				tri[k]);
		}
	}

	for (int layer_id = 0; layer_id < static_cast<int>(loops_by_layer.size()); ++layer_id) {
		if (layer_id >= static_cast<int>(flag_cut_layers_is_hole.size())
			|| flag_cut_layers_is_hole[layer_id] != -1) {
			continue;
		}

		const auto original_layer_it = follow_index.find(layer_id);
		if (original_layer_it == follow_index.end()) {
			std::cout << "[HybridManufacturing::CutMesh] Missing sorted cut-layer mapping for layer "
				<< layer_id << std::endl;
			continue;
		}
		const int original_layer_id = original_layer_it->second;
		const bool is_dependency_contact_face =
			layer_id < static_cast<int>(all_cut_layers_dependency_layer.size())
			&& all_cut_layers_dependency_layer[layer_id] > 0;

		for (auto outer_ring : loops_by_layer[layer_id]) {
			if (outer_ring.size() < 3) {
				continue;
			}
			if (CutRingSignedArea(outer_ring, all_slicer) < 0.0) {
				std::reverse(outer_ring.begin(), outer_ring.end());
			}

			std::vector<std::vector<int>> hole_rings;
			for (int hole_layer_id = 0; hole_layer_id < static_cast<int>(loops_by_layer.size()); ++hole_layer_id) {
				if (hole_layer_id >= static_cast<int>(flag_cut_layers_is_hole.size())
					|| flag_cut_layers_is_hole[hole_layer_id] != original_layer_id) {
					continue;
				}
				for (auto hole_ring : loops_by_layer[hole_layer_id]) {
					if (hole_ring.size() < 3
						|| !CutRingContainsPoint(
							outer_ring,
							all_slicer,
							all_slicer.positions[hole_ring.front()])) {
						continue;
					}
					if (CutRingSignedArea(hole_ring, all_slicer) > 0.0) {
						std::reverse(hole_ring.begin(), hole_ring.end());
					}
					hole_rings.push_back(std::move(hole_ring));
				}
			}

			std::vector<int> generated_cap_face_ids;
			if (contact_face_triangulation_mode == ContactFaceTriangulationMode::CGALRemesh) {
				vasco::contact_triangulation::TriangulateContactFacesWithCDT(
					outer_ring,
					hole_rings,
					all_slicer,
					generated_cap_face_ids,
					point_index_map);
			}
			else {
				using EarcutIndex = uint32_t;
				using EarcutPoint = std::array<double, 2>;
				std::vector<std::vector<EarcutPoint>> polygon;
				std::vector<int> flattened_vertex_ids;
				auto append_ring = [&](const std::vector<int>& ring) {
					polygon.emplace_back();
					polygon.back().reserve(ring.size());
					for (int vertex_id : ring) {
						const auto& point = all_slicer.positions[vertex_id];
						polygon.back().push_back({ point[0], point[1] });
						flattened_vertex_ids.push_back(vertex_id);
					}
				};
				append_ring(outer_ring);
				for (const auto& hole_ring : hole_rings) {
					append_ring(hole_ring);
				}

				const std::vector<EarcutIndex> indices = mapbox::earcut<EarcutIndex>(polygon);
				for (std::size_t i = 0; i + 2 < indices.size(); i += 3) {
					TRiangle triangle{
						flattened_vertex_ids[indices[i]],
						flattened_vertex_ids[indices[i + 1]],
						flattened_vertex_ids[indices[i + 2]]
					};
					const auto& a = all_slicer.positions[triangle[0]];
					const auto& b = all_slicer.positions[triangle[1]];
					const auto& c = all_slicer.positions[triangle[2]];
					if ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]) < 0.0) {
						std::swap(triangle[1], triangle[2]);
					}
					all_slicer.triangles.push_back(triangle);
					generated_cap_face_ids.push_back(static_cast<int>(all_slicer.triangles.size() - 1));
				}
			}

			all_cap_face_ids.insert(
				all_cap_face_ids.end(),
				generated_cap_face_ids.begin(),
				generated_cap_face_ids.end());
			if (is_dependency_contact_face) {
				id_contact_faces.insert(
					id_contact_faces.end(),
					generated_cap_face_ids.begin(),
					generated_cap_face_ids.end());
			}
		}
	}
}

void HybridManufacturing::RotateLayersForVisualization(
	vector<vector<cv::Point3d>>& all_layers,
	vector<vector<cv::Point3d>>& all_layers_contain,
	CutLayerHoles& all_layers_holes,
	const Eigen::Vector3d& vector_after,
	const Eigen::Vector3d& vector_before) const
{
	const Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vector_after, vector_before).toRotationMatrix();
	const auto rotate_layer = [&rotMatrix](CutLayer& layer) {
		for (cv::Point3d& point : layer) {
			const Eigen::Vector3d rotated = rotMatrix.inverse()
				* Eigen::Vector3d(point.x, point.y, point.z);
			point.x = rotated.x();
			point.y = rotated.y();
			point.z = rotated.z();
		}
	};
	for (int i = 0; i < all_layers.size(); i++) {
		rotate_layer(all_layers[i]);
		if (i < static_cast<int>(all_layers_contain.size())) {
			rotate_layer(all_layers_contain[i]);
		}
		if (i < static_cast<int>(all_layers_holes.size())) {
			for (CutLayer& hole : all_layers_holes[i]) {
				rotate_layer(hole);
			}
		}
	}
}

void HybridManufacturing::VisualizeCutLayers(
	const vector<vector<cv::Point3d>>& all_layers,
	const vector<vector<cv::Point3d>>& all_layers_contain,
	const CutLayerHoles& all_layers_holes,
	int height_of_beam_search,
	int cont_number_of_queue,
	int index_of_pre_node,
	bool judge_continue_additive,
	int id_continue,
	const Eigen::Vector3d& vector_after) const
{
	Visual Vis;
	cout << id_continue << endl;
	Vis.generateModelForRendering_5(all_layers, all_layers_contain, height_of_beam_search, cont_number_of_queue, file_name, index_of_pre_node, judge_continue_additive, id_continue);
	string vis_file(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "(" + to_string(index_of_pre_node) + ")_Layer");
	string vis_file_contain(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "(" + to_string(index_of_pre_node) + ")_Layer_contain");
	Eigen::Vector3f current_orientation = vector_after.cast<float>();
	if (open_vis_stair_case == true)
		visualize_layers_stair_case(vis_file, vis_file_contain, current_orientation, judge_continue_additive, id_continue);

	if (open_vis_stair_case == true && !map_S_and_vertex.empty()) {
		const std::string report_file = UnaccessibleVisPath(
			MakeAccessibilityDebugFileToken(vis_file) + "_subtractive_S_overlap.csv");
		std::ofstream report(report_file);
		if (report.is_open()) {
			report << std::setprecision(17);
			report << "layer_index,layer_z,overlap_s_count,s_id,vertex_id,x,y,z,z_diff,zero_ori_count,min_ori_count\n";

			const Eigen::Matrix3d rot_matrix =
				Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d(0, 0, 1), vector_after).toRotationMatrix();
			const Eigen::Matrix3d world_to_layer = rot_matrix.inverse();

			for (int layer_index = 0; layer_index < static_cast<int>(all_layers.size()); ++layer_index) {
				const auto& layer = all_layers[layer_index];
				if (layer.size() < 3) {
					continue;
				}

				std::vector<Point_2> layer_points;
				layer_points.reserve(layer.size());
				double layer_z = 0.0;
				for (const auto& point : layer) {
					const Eigen::Vector3d local_point = world_to_layer * Eigen::Vector3d(point.x, point.y, point.z);
					layer_points.emplace_back(local_point.x(), local_point.y());
					layer_z = local_point.z();
				}

				const double repair_eps = std::max(1e-9, dh * 1e-4);
				const double boundary_offset = std::max(1e-4, dh * 0.02);
				const LayerContainmentPolygon outer_polygon =
					PrepareLayerContainmentPolygon(layer_points, repair_eps);
				std::vector<LayerContainmentPolygon> hole_polygons;
				if (layer_index < static_cast<int>(all_layers_holes.size())) {
					hole_polygons.reserve(all_layers_holes[layer_index].size());
					for (const CutLayer& hole : all_layers_holes[layer_index]) {
						std::vector<Point_2> hole_points;
						hole_points.reserve(hole.size());
						for (const cv::Point3d& point : hole) {
							const Eigen::Vector3d local_point =
								world_to_layer * Eigen::Vector3d(point.x, point.y, point.z);
							hole_points.emplace_back(local_point.x(), local_point.y());
						}
						hole_polygons.push_back(
							PrepareLayerContainmentPolygon(hole_points, repair_eps));
					}
				}
				if (outer_polygon.points.size() < 3) {
					report << layer_index << "," << layer_z << ",invalid_polygon,-1,-1,0,0,0,0,0,0\n";
					continue;
				}

				std::vector<std::tuple<int, int, Eigen::Vector3d, double>> overlaps;
				for (const auto& entry : map_S_and_vertex) {
					const int s_id = entry.first;
					const int vertex_id = entry.second;
					if (vertex_id < 0 || vertex_id >= V.rows()) {
						continue;
					}

					const Eigen::Vector3d world_point(V(vertex_id, 0), V(vertex_id, 1), V(vertex_id, 2));
					const Eigen::Vector3d local_point = world_to_layer * world_point;
					const double z_diff = local_point.z() - layer_z;
					if (z_diff > 0.0
						&& z_diff <= dh
						&& IsPointInsidePreparedMaterialRegion(
							outer_polygon,
							hole_polygons,
							Point_2(local_point.x(), local_point.y()),
							boundary_offset)) {
						overlaps.emplace_back(s_id, vertex_id, world_point, z_diff);
					}
				}

				if (overlaps.empty()) {
					report << layer_index << "," << layer_z << ",0,-1,-1,0,0,0,0,0,0\n";
					continue;
				}

				for (const auto& overlap : overlaps) {
					const int s_id = std::get<0>(overlap);
					const int vertex_id = std::get<1>(overlap);
					const Eigen::Vector3d& point = std::get<2>(overlap);
					const double z_diff = std::get<3>(overlap);
					int zero_ori_count = 0;
					int min_ori_count = MAX_I;
					if (s_id >= 0 && s_id < static_cast<int>(ori_num_points_of_ori_in_all_the_area_S.size())) {
						for (int count : ori_num_points_of_ori_in_all_the_area_S[s_id]) {
							if (count == 0) {
								++zero_ori_count;
							}
							min_ori_count = std::min(min_ori_count, count);
						}
					}
					if (min_ori_count == MAX_I) {
						min_ori_count = -1;
					}
					report << layer_index << "," << layer_z << "," << overlaps.size() << ","
						<< s_id << "," << vertex_id << ","
						<< point.x() << "," << point.y() << "," << point.z() << ","
						<< z_diff << ","
						<< zero_ori_count << "," << min_ori_count << "\n";
				}
			}

			std::cout << "[VisualizeCutLayers] subtractive S overlap report saved to "
				<< report_file << std::endl;
		}
	}
}

void HybridManufacturing::Visualize_layer_segments(const std::vector<Layer>& layers) const
{
	vector<vector<cv::Point3d>> lines;
	for (const auto& layer : layers) {
		for (const auto& segment : layer.segments) {
			lines.emplace_back();
			auto& segment_line = lines.back();
			segment_line.emplace_back(segment.vertices[0].x, segment.vertices[0].y, layer.z);
			segment_line.emplace_back(segment.vertices[1].x, segment.vertices[1].y, layer.z);
		}
	}

	Visual vis;
	vis.generateModelForRendering_9(lines, file_name + "_katana_layers.obj");
}

void HybridManufacturing::Visualize_layer_polylines(const std::vector<Polylines>& layer_polylines) const
{
	vector<vector<cv::Point3d>> lines;
	for (const auto& layer : layer_polylines) {
		for (const auto& polyline : layer) {
			if (polyline.size() < 2) {
				continue;
			}
			lines.emplace_back();
			auto& polyline_line = lines.back();
			polyline_line.reserve(polyline.size());
			for (const auto& point : polyline) {
				polyline_line.emplace_back(point.x(), point.y(), point.z());
			}
		}
	}

	Visual vis;
	vis.generateModelForRendering_9(lines, file_name + "_cgal_polylines.obj");
}

Layer_Graph HybridManufacturing::BuildAdditiveLayerGraphWithSurfaceMesh(
	const Eigen::Vector3d& vector_after,
	int height_of_beam_search,
	int continue_node_id,
	const nozzle& the_nozzle,
	double& slicing_time,
	double& graph_time) const
{
	clock_t start_time_6 = clock();
	std::vector<SurfaceMeshSliceData> slices;
	BuildSurfaceMeshSlices(slices);

	clock_t end_time_6 = clock();
	slicing_time += double(end_time_6 - start_time_6) / CLOCKS_PER_SEC;

	std::vector<Polylines> contours;
	for (const auto& slice : slices) {
		contours.push_back(slice.contour_points);
		for (auto& contour : contours.back()) {
			contour.push_back(contour.front());
		}
	}
	//Visualize_layer_polylines(contours);
	std::vector<std::vector<std::vector<Eigen::Vector3d>> > surface_slice_points;
	std::vector<std::vector<std::vector<Eigen::Vector3d>> > surface_slice_points_contain;
	std::vector<std::vector<int>> surface_contour_parent_ids;
	surface_slice_points.reserve(slices.size());
	surface_slice_points_contain.reserve(slices.size());
	surface_contour_parent_ids.reserve(slices.size());
	for (const auto& slice : slices) {
		surface_slice_points.emplace_back();
		surface_slice_points_contain.emplace_back();
		surface_contour_parent_ids.push_back(slice.contour_parent_ids);
		for (const auto& contour : slice.contour_points) {
			std::vector<Eigen::Vector3d> contour_points;
			contour_points.reserve(contour.size());
			for (const auto& point : contour) {
				contour_points.push_back(Eigen::Vector3d(point.x(), point.y(), point.z()));
			}
			surface_slice_points.back().push_back(std::move(contour_points));
			surface_slice_points_contain.back().emplace_back();
		}
	}

	Data surface_data;
	surface_data.ReadData(
		surface_slice_points,
		surface_slice_points_contain,
		surface_contour_parent_ids);
	Layer_Graph layer_graph(surface_data);
	clock_t start_time_4 = clock();
	layer_graph.BuildLayerGraphFromSurfaceMeshSlices(slices, the_nozzle);
	PopulateLayerSelfSupportDiagnostics(
		slices,
		layer_graph);
	clock_t end_time_4 = clock();
	graph_time += double(end_time_4 - start_time_4) / CLOCKS_PER_SEC;

	return layer_graph;
}

std::vector<std::vector<int>> HybridManufacturing::EvaluateMergedPatchToolCollision(
	const Slicer_2& merged_patch,
	const std::vector<int>& merged_face_source_patch_id,
	cutter cutting_tool) const
{
	const int face_count = static_cast<int>(merged_patch.triangles.size());
	const int ori_count = static_cast<int>(sampling_subtractive.sample_points.size());

	// face x orientation: -1 表示无碰撞；>=0 表示碰撞面的最大 patch_index
	std::vector<std::vector<int>> min_collision_patch_matrix(face_count, std::vector<int>(ori_count, -1));
	if (face_count == 0 || ori_count == 0) {
		return min_collision_patch_matrix;
	}

	// 每个三角面的采样点（内心）
	std::vector<vasco::core::Vec3> sample_points(face_count);
	for (int i = 0; i < face_count; ++i) {
		const auto& tri = merged_patch.triangles[i];
		const auto& v1 = merged_patch.positions[tri[0]];
		const auto& v2 = merged_patch.positions[tri[1]];
		const auto& v3 = merged_patch.positions[tri[2]];

		double a = distance3d(v1, v2);
		double b = distance3d(v1, v3);
		double c = distance3d(v2, v3);

		vasco::core::Vec3 p{};
		if (a + b + c == 0.0) {
			p = v1;
		}
		else {
			for (int k = 0; k < 3; ++k) {
				p[k] = (a * v1[k] + b * v2[k] + c * v3[k]) / (a + b + c);
			}
		}
		sample_points[i] = p;
	}
	PrepareToolForCollision(cutting_tool);

	for (int ori = 0; ori < ori_count; ++ori) {
		std::vector<std::vector<Eigen::Vector3d>> temp_faces(face_count, std::vector<Eigen::Vector3d>(3));
		for (int i = 0; i < face_count; ++i) {
			for (int k = 0; k < 3; ++k) {
				int vid = merged_patch.triangles[i][k];
				temp_faces[i][k](0, 0) = merged_patch.positions[vid][0];
				temp_faces[i][k](1, 0) = merged_patch.positions[vid][1];
				temp_faces[i][k](2, 0) = merged_patch.positions[vid][2];
			}
		}

		std::vector<Eigen::Vector3d> temp_samples(face_count);
		for (int i = 0; i < face_count; ++i) {
			temp_samples[i](0, 0) = sample_points[i][0];
			temp_samples[i](1, 0) = sample_points[i][1];
			temp_samples[i](2, 0) = sample_points[i][2];
		}

		Eigen::Vector3d vectorBefore(0, 0, 1);
		Eigen::Vector3d vectorAfter(sampling_subtractive.sample_points[ori]);
		vectorAfter.normalize();

		Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vectorAfter).toRotationMatrix();

		for (int i = 0; i < face_count; ++i) {
			for (int k = 0; k < 3; ++k) {
				temp_faces[i][k] = rotMatrix.inverse() * temp_faces[i][k];
			}
			temp_samples[i] = rotMatrix.inverse() * temp_samples[i];
		}

		std::vector<double> max_z_of_faces(face_count, MIN_D);
		for (int i = 0; i < face_count; ++i) {
			for (int k = 0; k < 3; ++k) {
				max_z_of_faces[i] = std::max(max_z_of_faces[i], temp_faces[i][k](2, 0));
			}
		}

		std::vector<Eigen::Vector3d> normals(face_count);
		for (int i = 0; i < face_count; ++i) {
			const Eigen::Vector3d v1(temp_faces[i][0](0, 0), temp_faces[i][0](1, 0), temp_faces[i][0](2, 0));
			const Eigen::Vector3d v2(temp_faces[i][1](0, 0), temp_faces[i][1](1, 0), temp_faces[i][1](2, 0));
			const Eigen::Vector3d v3(temp_faces[i][2](0, 0), temp_faces[i][2](1, 0), temp_faces[i][2](2, 0));

			Eigen::Vector3d n = (v2 - v1).cross(v3 - v1);
			n.normalize();
			normals[i] = n;
		}

		// 计算：刀尖在第 i 面采样点、方向 ori 时，碰撞到的最大 patch_index
		for (int i = 0; i < face_count; ++i) {
			Eigen::Vector3d center_point;
			center_point.x() = temp_samples[i](0, 0) + cutting_tool.cylinder_r * normals[i].x();
			center_point.y() = temp_samples[i](1, 0) + cutting_tool.cylinder_r * normals[i].y();
			center_point.z() = temp_samples[i](2, 0) + cutting_tool.cylinder_r * normals[i].z();

			int max_patch_idx = std::numeric_limits<int>::min();
			bool has_collision = false;

			for (int ii = 0; ii < face_count; ++ii) {
				if (ii == i) {
					continue;
				}
				if (CheckToolCollisionWithCell(center_point, temp_faces[ii], max_z_of_faces[ii], cutting_tool, 30.0, 3.0)) {
					has_collision = true;
					if (ii < static_cast<int>(merged_face_source_patch_id.size())) {
						max_patch_idx = std::max(max_patch_idx, merged_face_source_patch_id[ii]);
					}
					else {
						max_patch_idx = std::max(max_patch_idx, -2);
					}
				}
			}

			// 无碰撞记为 -1；有碰撞则赋值为patch index + 1
			if (!has_collision) {
				min_collision_patch_matrix[i][ori] = -1;
			}
			else {
				min_collision_patch_matrix[i][ori] = max_patch_idx + 1;
			}
		}
	}

	//std::vector<int> min_collision_patch_per_face(face_count, -1);
	//for (int i = 0; i < face_count; ++i) {
	//	int row_min = std::numeric_limits<int>::max();
	//	for (int ori = 0; ori < ori_count; ++ori) {
	//		row_min = std::min(row_min, min_collision_patch_matrix[i][ori]);
	//	}
	//	if (row_min == std::numeric_limits<int>::max()) {
	//		row_min = -1;
	//	}
	//	min_collision_patch_per_face[i] = row_min;
	//}

	//// 可选输出
	//std::cout << "[Info] min_collision_patch_per_face size="
	//	<< min_collision_patch_per_face.size() << std::endl;

	//// 根据 min_collision_patch_per_face 输出带颜色的 OBJ（每个面一个颜色）
	//ExportMergedPatchFaceColorOBJ(
	//	merged_patch,
	//	min_collision_patch_per_face,
	//	".\\vis\\merged_patch_face_color_by_min_collision_patch.obj");

	return min_collision_patch_matrix;
}

void HybridManufacturing::ExportMergedPatchFaceColorOBJ(
	const Slicer_2& merged_patch,
	const std::vector<int>& max_collision_patch_per_face,
	const std::string& color_obj_file) const
{
	std::ofstream ofs(color_obj_file);
	if (!ofs.is_open()) {
		std::cout << "[Warn] cannot open file for writing: " << color_obj_file << std::endl;
		return;
	}

	int max_nonneg_label = -1;
	for (int v : max_collision_patch_per_face) {
		if (v >= 0) {
			max_nonneg_label = std::max(max_nonneg_label, v);
		}
	}
	std::cout << "[EvaluateMergedPatchToolCollision] max_nonneg_label = "
		<< max_nonneg_label << std::endl;

	auto color_from_label = [max_nonneg_label](int label) -> std::array<double, 3> {
		// -1: 无碰撞，固定颜色（浅灰）
		if (label < 0) {
			return { 0.80, 0.80, 0.80 };
		}

		// 非负: 越大越深（蓝色系）
		double t = 0.0;
		if (max_nonneg_label > 0) {
			t = static_cast<double>(label) / static_cast<double>(max_nonneg_label);
		}
		t = std::max(0.0, std::min(1.0, t));

		const std::array<double, 3> light = { 0.78, 0.88, 1.00 };
		const std::array<double, 3> dark = { 0.05, 0.20, 0.60 };

		std::array<double, 3> c;
		for (int k = 0; k < 3; ++k) {
			c[k] = light[k] * (1.0 - t) + dark[k] * t;
		}
		return c;
		};

	// 每个面写3个独立顶点，保证“面颜色”生效
	int v_count = 0;
	for (int i = 0; i < static_cast<int>(merged_patch.triangles.size()); ++i) {
		const auto& tri = merged_patch.triangles[i];
		const auto color = color_from_label(max_collision_patch_per_face[i]);

		for (int k = 0; k < 3; ++k) {
			const auto& p = merged_patch.positions[tri[k]];
			ofs << "v " << p[0] << " " << p[1] << " " << p[2] << " "
				<< color[0] << " " << color[1] << " " << color[2] << "\n";
		}

		ofs << "f " << (v_count + 1) << " " << (v_count + 2) << " " << (v_count + 3) << "\n";
		v_count += 3;
	}

	std::cout << "[Info] wrote colored OBJ: " << color_obj_file << std::endl;
}

Slicer_2 HybridManufacturing::MergeBlockPatchesWithDedup(
	const std::vector<std::string>& patch_files,
	std::vector<int>& merged_vertex_source_patch_id,
	std::vector<int>& merged_face_source_patch_id,
	double merge_eps) const
{
	Slicer_2 merged;
	merged_vertex_source_patch_id.clear();
	merged_face_source_patch_id.clear();

	struct QuantizedVertexKey {
		long long x;
		long long y;
		long long z;
		bool operator<(const QuantizedVertexKey& rhs) const {
			if (x != rhs.x) return x < rhs.x;
			if (y != rhs.y) return y < rhs.y;
			return z < rhs.z;
		}
	};

	if (merge_eps <= 0.0) {
		merge_eps = 1e-6;
	}
	const double inv_eps = 1.0 / merge_eps;
	auto make_vertex_key = [inv_eps](const vasco::Slicer::Vec3& p) {
		return QuantizedVertexKey{
			static_cast<long long>(std::llround(p[0] * inv_eps)),
			static_cast<long long>(std::llround(p[1] * inv_eps)),
			static_cast<long long>(std::llround(p[2] * inv_eps))
		};
		};

	std::map<QuantizedVertexKey, int> vertex_map;
	std::map<std::array<int, 3>, int> face_key_to_merged_face_index;

	for (int patch_list_index = 0; patch_list_index < static_cast<int>(patch_files.size()); ++patch_list_index) {
		const int patch_index = patch_list_index + 1;
		const std::string& patch_file = patch_files[patch_list_index];

		Slicer_2 patch;
		if (!patch.load(patch_file)) {
			std::cout << "[Warn] cannot load ancestor patch file: " << patch_file << std::endl;
			continue;
		}
		const std::size_t original_vertex_count = patch.positions.size();
		std::size_t removed_vertex_count = 0;
		if (!CompactPatchToReferencedVertices(patch, removed_vertex_count)) {
			std::cout << "[Warn] skip ancestor patch with invalid face indices: "
				<< patch_file << std::endl;
			continue;
		}

		std::cout << "[Info] merge ancestor patch " << patch_index
			<< ": " << patch_file
			<< ", V=" << original_vertex_count
			<< ", referenced_V=" << patch.positions.size()
			<< ", removed_unreferenced_V=" << removed_vertex_count
			<< ", F=" << patch.triangles.size() << std::endl;

		std::vector<int> local_to_global(patch.positions.size(), -1);

		for (int i = 0; i < static_cast<int>(patch.positions.size()); ++i) {
			const auto key = make_vertex_key(patch.positions[i]);
			auto it = vertex_map.find(key);

			if (it == vertex_map.end()) {
				const int new_index = static_cast<int>(merged.positions.size());
				merged.positions.push_back(patch.positions[i]);
				vertex_map.insert({ key, new_index });
				local_to_global[i] = new_index;
				merged_vertex_source_patch_id.push_back(patch_index);
			}
			else {
				local_to_global[i] = it->second;
			}
		}

		for (const auto& tri : patch.triangles) {
			int a = local_to_global[tri[0]];
			int b = local_to_global[tri[1]];
			int c = local_to_global[tri[2]];

			if (a == b || b == c || a == c) {
				continue;
			}

			std::array<int, 3> key = { a, b, c };
			std::sort(key.begin(), key.end());

			auto fit = face_key_to_merged_face_index.find(key);
			if (fit == face_key_to_merged_face_index.end()) {
				const int new_face_index = static_cast<int>(merged.triangles.size());
				merged.triangles.push_back({ a, b, c });
				face_key_to_merged_face_index.insert({ key, new_face_index });
				merged_face_source_patch_id.push_back(patch_index);
			}
		}
	}

	return merged;
}

Slicer_2 HybridManufacturing::MergeBlockPatchesWithDedup(
	int max_patch_index,
	std::vector<int>& merged_vertex_source_patch_id,
	std::vector<int>& merged_face_source_patch_id,
	double merge_eps) const
{
	Slicer_2 merged;
	merged_vertex_source_patch_id.clear();
	merged_face_source_patch_id.clear();

	struct QuantizedVertexKey {
		long long x;
		long long y;
		long long z;
		bool operator<(const QuantizedVertexKey& rhs) const {
			if (x != rhs.x) return x < rhs.x;
			if (y != rhs.y) return y < rhs.y;
			return z < rhs.z;
		}
	};

	const double inv_eps = 1.0 / merge_eps;
	auto make_vertex_key = [inv_eps](const vasco::Slicer::Vec3& p) {
		return QuantizedVertexKey{
			static_cast<long long>(std::llround(p[0] * inv_eps)),
			static_cast<long long>(std::llround(p[1] * inv_eps)),
			static_cast<long long>(std::llround(p[2] * inv_eps))
		};
		};

	std::map<QuantizedVertexKey, int> vertex_map;
	std::map<std::array<int, 3>, int> face_key_to_merged_face_index;

	for (int patch_index = 1; patch_index <= max_patch_index; ++patch_index) {
		Slicer_2 patch;
		const std::string patch_file = PatchVisPath(
			"block_patch-" + to_string(patch_index) + "_0.obj");

		if (!patch.load(patch_file)) {
			std::cout << "[Warn] cannot load patch file: " << patch_file << std::endl;
			continue;
		}
		const std::size_t original_vertex_count = patch.positions.size();
		std::size_t removed_vertex_count = 0;
		if (!CompactPatchToReferencedVertices(patch, removed_vertex_count)) {
			std::cout << "[Warn] skip patch with invalid face indices: "
				<< patch_file << std::endl;
			continue;
		}
		std::cout << "[Info] compact patch " << patch_index
			<< ": " << patch_file
			<< ", V=" << original_vertex_count
			<< ", referenced_V=" << patch.positions.size()
			<< ", removed_unreferenced_V=" << removed_vertex_count
			<< ", F=" << patch.triangles.size() << std::endl;

		std::vector<int> local_to_global(patch.positions.size(), -1);

		for (int i = 0; i < static_cast<int>(patch.positions.size()); ++i) {
			const auto key = make_vertex_key(patch.positions[i]);
			auto it = vertex_map.find(key);

			if (it == vertex_map.end()) {
				const int new_index = static_cast<int>(merged.positions.size());
				merged.positions.push_back(patch.positions[i]);
				vertex_map.insert({ key, new_index });
				local_to_global[i] = new_index;

				// 单归属：首次出现即最小 patch_index
				merged_vertex_source_patch_id.push_back(patch_index);
			}
			else {
				// 已存在：保持原归属（更小 patch_index）
				local_to_global[i] = it->second;
			}
		}

		for (const auto& tri : patch.triangles) {
			int a = local_to_global[tri[0]];
			int b = local_to_global[tri[1]];
			int c = local_to_global[tri[2]];

			if (a == b || b == c || a == c) {
				continue;
			}

			std::array<int, 3> key = { a, b, c };
			std::sort(key.begin(), key.end());

			auto fit = face_key_to_merged_face_index.find(key);
			if (fit == face_key_to_merged_face_index.end()) {
				const int new_face_index = static_cast<int>(merged.triangles.size());
				merged.triangles.push_back({ a, b, c });
				face_key_to_merged_face_index.insert({ key, new_face_index });

				// 单归属：首次出现即最小 patch_index
				merged_face_source_patch_id.push_back(patch_index);
			}
			// 已存在：保持原归属（更小 patch_index）
		}
	}

	return merged;
}


void HybridManufacturing::CutMesh(
	CutLayerVector all_layers,
	CutLayerVector all_layers_contain,
	CutLayerHoles all_layers_holes,
	CutLayerVector all_cut_layers,
	Eigen::Vector3d vector_after,
	int height_of_beam_search,
	int cont_number_of_queue,
	int index_of_pre_node,
	vector<int> all_cut_layers_dependency_layer,
	bool& jud_outer_beam_search_terminate,
	vector<TRiangle>& current_remove_triangles,
	Slicer_2& current_slicer,
	bool judge_continue_additive,
	bool flag_is_continue_block,
	int pre_cont_number_of_queue,
	bool& jud_error,
	int id_node,
	int id_continue,
	vector<int> flag_cut_layers_is_hole,
	const vector<bool>& judge_S_be_searched)
{
	std::cout << "[HybridManufacturing::CutMesh] height_of_beam_search: " << height_of_beam_search
		<< " cont_number_of_queue: " << cont_number_of_queue
		<< " index_of_pre_node: " << index_of_pre_node
		<< " id_node: " << id_node
		<< " id_continue: " << id_continue
		<< std::endl;
	bool using_solid_model = true;

	std::vector<bool> jud_triangle_have_been_added;	//三角形是否被添加进remove_triangles
	vector<vector<int>> cutting_plane_points;	//每个切割平面对应的切割点
	cutting_plane_points.resize(all_cut_layers.size());
	vector<vector<pair<int, int>>> cutting_plane_edges;	//每个切割平面对应的切割边 <顶点index,顶点index>
	cutting_plane_edges.resize(all_cut_layers.size());

	if (all_cut_layers.size() == 0) {
		std::cout << "[HybridManufacturing::CutMesh]No cut layers!" << std::endl;
		return;
	}
	//////sort cut layers//////
	map<int, int> follow_index;	//记录排序前后cut layers的对应关系 <排序后index,排序前index>
	SortCutLayersByHeight(all_cut_layers, all_cut_layers_dependency_layer, flag_cut_layers_is_hole, follow_index);
	///////////////////////////
	//cout << "a" << id_continue;
	Slicer_2 slicer = LoadSlicerForCutMesh(flag_is_continue_block, height_of_beam_search, index_of_pre_node, id_continue);


	if (flag_is_continue_block == true) {
		height_of_beam_search--;
		cont_number_of_queue = pre_cont_number_of_queue;
	}

	//need rotate first//
	Eigen::Vector3d vectorBefore(0, 0, 1);
	RotateSlicerPositions(slicer, vectorBefore, vector_after);
	/////////////////////

	//layer visualization//
	RotateLayersForVisualization(
		all_layers,
		all_layers_contain,
		all_layers_holes,
		vector_after,
		vectorBefore);
	cout << "b" << endl;
	VisualizeCutLayers(
		all_layers,
		all_layers_contain,
		all_layers_holes,
		height_of_beam_search,
		cont_number_of_queue,
		index_of_pre_node,
		judge_continue_additive,
		id_continue,
		vector_after);
	////////////////////////
	cout << "B" << endl;
	//--------------------cut-----------------------//
	slicer.normal[0] = 0;
	slicer.normal[1] = 0;
	slicer.normal[2] = 1;

	slicer.jud_plane.resize(slicer.triangles.size());
	for (int i = 0; i < slicer.triangles.size(); i++)
		slicer.jud_plane[i] = false;
	clock_t start_time_3, end_time_3;

	for (int i = 0; i < all_cut_layers.size(); i++) {
		slicer.origin[0] = 0;
		slicer.origin[1] = 0;
		slicer.origin[2] = all_cut_layers[i][0].z;
		slicer.cut();
		//cout << "ok";
		cout << "i " << i << " slicer.origin[2] " << slicer.origin[2] << endl;
	}

	std::vector<TRiangle> ori_triangle = slicer.triangles;	//保存切割后的所有三角形
	jud_triangle_have_been_added.resize(slicer.triangles.size());	//为ori_triangle中的每个三角形分配一个标记，表示是否被添加进remove_triangles
	for (int i = 0; i < slicer.triangles.size(); i++)
		jud_triangle_have_been_added[i] = false;

	Slicer_2 all_slicer;	//all_slicer保存切割后的所有三角形和顶点

	all_slicer.positions = slicer.positions;
	all_slicer.triangles = ori_triangle;


	all_slicer.save(file_name + "-after_cut-" + std::to_string(height_of_beam_search) + "EEEERQ.obj");

#if 0 // Legacy height-sorted, vertex-adjacent removal propagation.
	//--------------------save candidate_triangles-----------------------//
	int current_index = 0;	//
	std::vector<VEctor> min_z_point;	//每个候选三角形中z值最小的顶点坐标
	std::vector<double> min_z_triangle;	//每个候选三角形中z值最小的顶点的z值
	vector<int> index_of_min_point_in_triangle;	//每个候选三角形中z值最小的顶点在三角形中的索引
	std::vector<TRiangle> candidate_triangles;	//候选三角形集合
	std::vector<int> id_candidate_triangles;	//候选三角形在all_slicer.triangles中的索引
	std::vector<int> id_triangles;	//候选三角形在slicer.triangles中的索引
	std::vector<TRiangle> remove_triangles;	//最终被移除的三角形集合
	std::vector<int> id_remove_triangles;	//最终被移除的三角形在all_slicer.triangles中的索引
	candidate_triangles.clear();
	id_triangles.clear();
	double min_z_all_cut_layers = 9999999;	//所有切割层中z值最小的值
	for (int t = 0; t < all_cut_layers.size(); t++)
		min_z_all_cut_layers = min(min_z_all_cut_layers, all_cut_layers[t][0].z);
	start_time_3 = clock();
	for (int i = 0; i < all_slicer.triangles.size(); i++) {	//遍历所有切割后的三角形
		double temp_min_z_triangle = 9999999;	//当前三角形中z值最小的顶点的z值
		for (int k = 0; k < 3; k++) {
			temp_min_z_triangle = min(all_slicer.positions[all_slicer.triangles[i][k]][2], temp_min_z_triangle);
		}
		if (temp_min_z_triangle + 0.001 >= min_z_all_cut_layers) {	//如果当前三角形中z值最小的顶点的z值大于等于所有切割层中z值最小的值，则将该三角形加入候选三角形集合
			candidate_triangles.push_back(all_slicer.triangles[i]);	//加入候选三角形集合
			id_candidate_triangles.push_back(i);	//记录该候选三角形在all_slicer.triangles中的索引
			for (int k = i; k < slicer.triangles.size(); k++)
				if (slicer.triangles[k] == all_slicer.triangles[i]) {	//找到该候选三角形在slicer.triangles中的索引
					id_triangles.push_back(k);
					break;
				}
		}
	}
	end_time_3 = clock();
	for (int i = 0; i < candidate_triangles.size(); i++) {
		double temp_min_z_triangle = 9999999;	//当前候选三角形中z值最小的顶点的z值
		VEctor temp_min_z_point = all_slicer.positions[candidate_triangles[i][0]];
		int temp_index_of_min_point_in_triangle;
		for (int k = 0; k < 3; k++) {
			if (all_slicer.positions[candidate_triangles[i][k]][2] < temp_min_z_triangle) {	//找到当前候选三角形中z值最小的顶点
				temp_min_z_triangle = all_slicer.positions[candidate_triangles[i][k]][2];	//更新z值最小的顶点的z值
				temp_min_z_point = all_slicer.positions[candidate_triangles[i][k]];	//更新z值最小的顶点的坐标
				temp_index_of_min_point_in_triangle = k;	//更新z值最小的顶点在三角形中的索引
			}
		}
		min_z_triangle.push_back(temp_min_z_triangle);	//将当前候选三角形中z值最小的顶点的z值加入min_z_triangle
		min_z_point.push_back(temp_min_z_point);	//将当前候选三角形中z值最小的顶点的坐标加入min_z_point
		index_of_min_point_in_triangle.push_back(temp_index_of_min_point_in_triangle);	//将当前候选三角形中z值最小的顶点在三角形中的索引加入index_of_min_point_in_triangle
	}

	cout << "()()()(" << double(end_time_3 - start_time_3) / CLOCKS_PER_SEC << endl;
	clock_t start_time_2, end_time_2;
	start_time_2 = clock();
	for (int i = 0; i < candidate_triangles.size(); i++) {	//对候选三角形按照min_z_triangle进行排序
		for (int j = i + 1; j < candidate_triangles.size(); j++) {
			if (min_z_triangle[i] > min_z_triangle[j]) {
				swap(candidate_triangles[i], candidate_triangles[j]);	//交换候选三角形
				swap(id_candidate_triangles[i], id_candidate_triangles[j]);	//交换候选三角形在all_slicer.triangles中的索引
				swap(id_triangles[i], id_triangles[j]);	//交换候选三角形在slicer.triangles中的索引
				swap(min_z_triangle[i], min_z_triangle[j]);	//交换min_z_triangle
				swap(min_z_point[i], min_z_point[j]);	//交换min_z_point
				swap(index_of_min_point_in_triangle[i], index_of_min_point_in_triangle[j]);	//交换index_of_min_point_in_triangle
			}
		}
	}
	end_time_2 = clock();
	cout << "()()()(" << double(end_time_2 - start_time_2) / CLOCKS_PER_SEC << endl;
	//slicer.triangles = candidate_triangles;

	//--------------------save OPP_triangles one by one-----------------------//
	vector<int> save_current_index;
	for (int t = 0; t < all_cut_layers.size(); t++) {	//枚举all_cut_layers中的每一层
		current_index = 0;   //该方式也许比较慢
		double boundary_bottom = 999999, boundary_left = 999999, boundary_top = -999999, boundary_right = -999999;	//当前切割层的边界
		//Point_2* points = new Point_2[all_cut_layers[t].size()];
		vector<Point_2> points;
		points.resize(all_cut_layers[t].size());
		for (int i = 0; i < all_cut_layers[t].size(); i++) {	//枚举all_cut_layers[t]中的点，计算当前切割层的二维AABB
			boundary_top = std::max(boundary_top, all_cut_layers[t][i].y);
			boundary_bottom = std::min(boundary_bottom, all_cut_layers[t][i].y);
			boundary_right = std::max(boundary_right, all_cut_layers[t][i].x);
			boundary_left = std::min(boundary_left, all_cut_layers[t][i].x);
			Point_2 temp_point(all_cut_layers[t][i].x, all_cut_layers[t][i].y);
			points[i] = temp_point;
		}
		const LayerContainmentPolygon prepared_cut_layer_polygon =
			PrepareLayerContainmentPolygon(points, 1e-8);
		cv::Point2d current_triangle_point;
		cv::Point2d current_layer_point;
		for (; current_index < candidate_triangles.size(); current_index++) {	//枚举所有候选三角形
			if (abs(min_z_triangle[current_index] - all_cut_layers[t][0].z) > 0.0001) {
				if (min_z_triangle[current_index] > all_cut_layers[t][0].z) {	//如果当前候选三角形中z值最小的顶点的z值大于当前切割层的z值，则跳出循环
					break;
				}
			}
			else {	//如果当前候选三角形中最小z值等于当前切割层的z值,才进行是否加入remove_triangles的判断
				CutEdgeKey cut_boundary_edge;
				const bool jud_is_boundary_point =
					TryGetUpperCutBoundaryEdge(
						all_slicer,
						candidate_triangles[current_index],
						all_cut_layers[t][0].z,
						cut_boundary_edge)
					&& CutBoundaryEdgeBelongsToLayer(
						all_slicer,
						cut_boundary_edge,
						prepared_cut_layer_polygon,
						1e-3);
				/*if (min_z_point[current_index][0] + 0.01 >= boundary_left && min_z_point[current_index][0] - 0.01 <= boundary_right &&
					min_z_point[current_index][1] + 0.01 >= boundary_bottom && min_z_point[current_index][1] - 0.01 <= boundary_top) {
					jud_is_boundary_point = true;
					break;
				}*/
				//}

				if (jud_is_boundary_point == true) {
					auto layer_begin = points.begin();
					auto layer_end = points.end();
					//if (check_inside_2(Point_2(current_triangle_point.x, current_triangle_point.y), points, points + all_cut_layers[t].size(), K())) {
					current_triangle_point.x = min_z_point[current_index][0];
					current_triangle_point.y = min_z_point[current_index][1];
					/*if (check_inside_2(Point_2(current_triangle_point.x, current_triangle_point.y), layer_begin, layer_end, K())) {*/

					if (1) {
						int cont_cutting_points = 0;
						int temp_left = 0, temp_right = 0;
						for (int k = 0; k < 3; k++) {
							if (abs(all_cut_layers[t][0].z - all_slicer.positions[candidate_triangles[current_index][k]][2]) < 0.0001) {	//判断当前顶点是否在切割平面上（前面已经保证了z值最小的顶点不大于切割平面z值）
								cont_cutting_points++;
								if (cont_cutting_points == 1)
									temp_left = candidate_triangles[current_index][k];
								else if (cont_cutting_points == 2) {
									temp_right = candidate_triangles[current_index][k];
									cutting_plane_edges[t].push_back(make_pair(temp_left, temp_right));
									break;
								}
								//cutting_plane_points[t].push_back(candidate_triangles[current_index][index_of_min_point_in_triangle[current_index]]);
							}
						}
						remove_triangles.push_back(candidate_triangles[current_index]);	//将当前候选三角形加入remove_triangles
						id_remove_triangles.push_back(id_candidate_triangles[current_index]);	//将当前候选三角形在all_slicer.triangles中的索引加入id_remove_triangles
						/*if (height_of_beam_search == 2 && id_triangles[current_index] == 10297)
							cout << "********************************** " << all_cut_layers[t][0].z << endl;*/

						save_current_index.push_back(current_index);	//记录当前候选三角形在candidate_triangles中的索引到save_current_index
						jud_triangle_have_been_added[id_triangles[current_index]] = true;	//jud_triangle_have_been_added标记当前候选三角形已被添加进remove_triangles
					}
					else {
						for (int j = 0; j < all_cut_layers[t].size(); j++) {
							current_layer_point.x = all_cut_layers[t][j].x;	//
							current_layer_point.y = all_cut_layers[t][j].y;
							current_triangle_point.x = min_z_point[current_index][0];
							current_triangle_point.y = min_z_point[current_index][1];
							//if ((jud_is_boundary_point == false && distance_2d(current_layer_point, current_triangle_point) < 0.002) || (jud_is_boundary_point == true)) {   //0.002
							if (distance2d(current_layer_point, current_triangle_point) < 0.1) {  //4.0
								int cont_cutting_points = 0;
								int temp_left = 0, temp_right = 0;
								for (int k = 0; k < 3; k++) {
									if (abs(all_cut_layers[t][0].z - all_slicer.positions[candidate_triangles[current_index][k]][2]) < 0.0001) {
										cont_cutting_points++;
										if (cont_cutting_points == 1)
											temp_left = candidate_triangles[current_index][k];
										else if (cont_cutting_points == 2) {
											temp_right = candidate_triangles[current_index][k];
											cutting_plane_edges[t].push_back(make_pair(temp_left, temp_right));
											break;
										}
										//cutting_plane_points[t].push_back(candidate_triangles[current_index][index_of_min_point_in_triangle[current_index]]);
									}
								}
								remove_triangles.push_back(candidate_triangles[current_index]);
								id_remove_triangles.push_back(id_candidate_triangles[current_index]);
								/*if (height_of_beam_search == 2 && id_triangles[current_index] == 10297)
									cout << "********************************** " << all_cut_layers[t][0].z << endl;*/

								save_current_index.push_back(current_index);
								jud_triangle_have_been_added[id_triangles[current_index]] = true;
								break;
							}
						}
					}
				}
			}

		}
	}
	current_index = 0;

	while (1) {
		bool flag_break = false;
		for (int i = 0; i < remove_triangles.size(); i++)
		{
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					if (jud_triangle_have_been_added[id_triangles[current_index]] == false && candidate_triangles[current_index][j] == remove_triangles[i][k]
						&& min_z_triangle[current_index] >= min_z_triangle[save_current_index[i]]) {	//如果当前候选三角形的某个顶点与remove_triangles[i]的某个顶点相同，且当前候选三角形中z值最小的顶点的z值大于等于remove_triangles[i]中z值最小的顶点的z值
						//if (height_of_beam_search != 2)
						remove_triangles.push_back(candidate_triangles[current_index]);	//将当前候选三角形加入remove_triangles
						id_remove_triangles.push_back(id_candidate_triangles[current_index]);
						save_current_index.push_back(current_index);
						jud_triangle_have_been_added[id_triangles[current_index]] = true;
						flag_break = true;
						break;
					}
				}
				if (flag_break == true)
					break;
			}
			if (flag_break == true)
				break;
		}
		current_index++;
		if (current_index >= candidate_triangles.size())
			break;
	}

#endif

	const std::vector<int> seed_remove_face_ids = CollectCutRemovalSeedFaces(
		all_slicer,
		all_cut_layers,
		cutting_plane_edges);
	jud_triangle_have_been_added = ExpandRemovedFacesAcrossNonCutEdges(
		all_slicer,
		seed_remove_face_ids,
		cutting_plane_edges);

	std::vector<TRiangle> remove_triangles;
	std::vector<int> id_remove_triangles;
	remove_triangles.reserve(all_slicer.triangles.size());
	id_remove_triangles.reserve(all_slicer.triangles.size());
	for (int face_id = 0; face_id < static_cast<int>(all_slicer.triangles.size()); ++face_id) {
		if (!jud_triangle_have_been_added[face_id]) {
			continue;
		}
		remove_triangles.push_back(all_slicer.triangles[face_id]);
		id_remove_triangles.push_back(face_id);
	}
	std::cout << "[HybridManufacturing::CutMesh] Edge-connected removal selection: "
		<< "seed_faces=" << seed_remove_face_ids.size()
		<< ", removed_faces=" << remove_triangles.size()
		<< std::endl;
	if (remove_triangles.size() >= all_slicer.triangles.size()) {
		std::cerr << "[HybridManufacturing::CutMesh] Reject cut because removal selection consumed the entire mesh: "
			<< "node=" << id_node
			<< ", removed_faces=" << remove_triangles.size()
			<< ", total_faces=" << all_slicer.triangles.size()
			<< std::endl;
		jud_error = true;
		current_remove_triangles.clear();
		current_slicer.clear();
		return;
	}

	vasco::Slicer temp_slicer_1;
	temp_slicer_1.clear();
	temp_slicer_1.positions = all_slicer.positions;
	for (int i = 0; i < id_remove_triangles.size(); i++) {
		temp_slicer_1.triangles.push_back(all_slicer.triangles[id_remove_triangles[i]]);
	}
	//cout << "remove faces " << temp_slicer_1.triangles.size() << endl;

	//RotateSlicerPositions(temp_slicer_1, vector_after, vectorBefore);
	temp_slicer_1.save(file_name + "-all2-" + std::to_string(height_of_beam_search) + "QWERQ.obj");

	std::vector<std::vector<Eigen::Vector3d>> all_cut_layers_v3d;
	for (int i = 0; i < all_cut_layers.size(); i++) {
		std::vector<Eigen::Vector3d> temp_layer;
		for (int j = 0; j < all_cut_layers[i].size(); j++) {
			Eigen::Vector3d temp_point(all_cut_layers[i][j].x, all_cut_layers[i][j].y, all_cut_layers[i][0].z);
			temp_layer.push_back(temp_point);
		}
		all_cut_layers_v3d.push_back(temp_layer);
	}
	vasco::io::outputContourToObj(all_cut_layers_v3d, Eigen::Vector3d(0.5, 0.5, 0.0), file_name + "-all_cut_layers-" + std::to_string(height_of_beam_search) + "QWERQ.obj");


	// Legacy per-triangle residual-face removal stays disabled because it creates off-plane holes.
	// for (int i = 0; i < jud_triangle_have_been_added.size(); i++) {
	// 	if (jud_triangle_have_been_added[i] == false) {
	// 		int id_layer;
	// 		double min_dis = 99999999;
	// 		for (int t = 0; t < all_cut_layers.size(); t++) {
	// 			if (all_cut_layers[t][0].z - slicer.positions[slicer.triangles[i][0]][2] <= 2 * dh && all_cut_layers[t][0].z - slicer.positions[slicer.triangles[i][0]][2] >= -0.001) {	//如果当前切割层的z值与当前三角形中z值最小的顶点的z值之差在[-0.001,2*dh]范围内
	// 				for (int j = 0; j < all_cut_layers[t].size(); j++) {
	// 					cv::Point3d current_triangle_point(slicer.positions[slicer.triangles[i][0]][0], slicer.positions[slicer.triangles[i][0]][1], slicer.positions[slicer.triangles[i][0]][2]);
	// 					cv::Point3d current_layer_point(all_cut_layers[t][j].x, all_cut_layers[t][j].y, all_cut_layers[t][0].z);

	// 					double distance = distance3d(current_triangle_point, current_layer_point);
	// 					if (distance < min_dis) {
	// 						min_dis = distance;
	// 						id_layer = t;
	// 					}
	// 				}
	// 			}
	// 		}
	// 		if (min_dis < Dis && all_cut_layers_dependency_layer[id_layer] == 0) {
	// 			cout << "还真的能到这个地方啊 - 残余面片 min_dis " << min_dis << "Dis " << Dis << std::endl;
	// 			remove_triangles.push_back(slicer.triangles[i]);
	// 			id_remove_triangles.push_back(i);
	// 			//jud_triangle_have_been_added[i] = true;
	// 		}

	// 	}
	// }




	/////////////////////////删除残余面片////////////////////////////
	//建立面片邻接关系
#if 0 // Replaced by linear-time edge union below.
	vector<vector<int>> adjacent_faces(slicer.triangles.size());
	//建立不在remove_triangles中的面片之间的邻接关系
	for (int i = 0; i < slicer.triangles.size(); i++) {
		if (jud_triangle_have_been_added[i] == false) {
			for (int j = i + 1; j < slicer.triangles.size(); j++) {
				if (jud_triangle_have_been_added[j] == false) {
					int cont_same_point = 0;
					for (int k = 0; k < 3; k++) {
						for (int l = 0; l < 3; l++) {
							if (slicer.triangles[i][k] == slicer.triangles[j][l]) {
								cont_same_point++;
								break;
							}
						}
					}
					if (cont_same_point >= 2) {	//如果两个面片有两个以上的公共顶点，则认为它们是邻接的
						adjacent_faces[i].push_back(j);
						adjacent_faces[j].push_back(i);
					}
				}
			}
		}
	}
	//分区
	vector<bool> visited(slicer.triangles.size(), false);
	vector<vector<int>> save_faces(0);
	for (int i = 0; i < slicer.triangles.size(); ++i)
	{
		if (jud_triangle_have_been_added[i] || visited[i]) {
			continue; // 已标记或已访问，跳过
		}

		visited[i] = true;

		// 当前连通分量的索引集合
		std::vector<int> component;
		component.push_back(i);

		// 标准队列进行 BFS
		std::queue<int> q;
		q.push(i);
		while (!q.empty())
		{
			const int u = q.front();
			q.pop();

			// 遍历 u 的邻接面
			const auto& adj = adjacent_faces[u];
			for (int v : adj)
			{
				if (!visited[v] && !jud_triangle_have_been_added[v])
				{
					visited[v] = true;
					q.push(v);
					component.push_back(v);
				}
			}
		}

		// 将本连通分量保存（对应之前的 save_faces.push_back(...)）
		save_faces.push_back(std::move(component));
	}

	// 移除过小的残余面片分量
	for (const auto& comp : save_faces)
	{
		cout << "&&*&& " << comp.size() << endl;
		if (comp.size() < 80)
		{
			for (int idx : comp)
			{
				remove_triangles.push_back(slicer.triangles[idx]);
				id_remove_triangles.push_back(idx);
				jud_triangle_have_been_added[idx] = true;
			}
		}
	}
#endif
	const std::vector<std::vector<int>> remaining_face_components =
		CollectUnselectedFaceComponents(slicer, jud_triangle_have_been_added);
	for (const auto& component : remaining_face_components) {
		std::cout << "&&*&& " << component.size() << std::endl;
		if (component.size() >= 80) {
			continue;
		}
		for (int face_id : component) {
			if (jud_triangle_have_been_added[face_id]) {
				continue;
			}
			remove_triangles.push_back(slicer.triangles[face_id]);
			id_remove_triangles.push_back(face_id);
			jud_triangle_have_been_added[face_id] = true;
		}
	}
	if (remove_triangles.size() >= all_slicer.triangles.size()) {
		std::cerr << "[HybridManufacturing::CutMesh] Reject cut because removal plus residual cleanup consumed the entire mesh: "
			<< "node=" << id_node
			<< ", removed_faces=" << remove_triangles.size()
			<< ", total_faces=" << all_slicer.triangles.size()
			<< std::endl;
		jud_error = true;
		current_remove_triangles.clear();
		current_slicer.clear();
		return;
	}
	//////////////////////////////////////////////////////////////////



	cout << "c" << endl;
	//current_remove_triangles = remove_triangles;

	current_remove_triangles = remove_triangles; // current_remove_triangles先不进行筛选，直接等于remove_triangles，后面再剔除非表面的面片
	cout << "ccc" << endl;

	//将slicer旋转回原始位置
	RotateSlicerPositions(slicer, vector_after, vectorBefore);
	current_slicer = slicer;

	size_t face_cnt = Normals.rows();
	std::cout << "face_cnt " << face_cnt << std::endl;
	current_remove_triangles = FilterSurfaceRemoveTriangles(slicer, remove_triangles);
	const auto unresolved_area_s_hits = FindUnresolvedAreaSOnRemovedSurface(
		current_slicer,
		current_remove_triangles,
		V,
		map_S_and_vertex,
		judge_S_be_searched);
	if (!unresolved_area_s_hits.empty()) {
		std::cerr << "[HybridManufacturing::CutMesh] Reject cut because the removed surface contains unresolved subtractive-inaccessible vertices: "
			<< "node=" << id_node
			<< ", hit_count=" << unresolved_area_s_hits.size()
			<< std::endl;
		const std::size_t report_count = std::min<std::size_t>(unresolved_area_s_hits.size(), 20);
		for (std::size_t hit_index = 0; hit_index < report_count; ++hit_index) {
			const int s_id = unresolved_area_s_hits[hit_index].first;
			const int vertex_id = unresolved_area_s_hits[hit_index].second;
			std::cerr << "  area_S=" << s_id
				<< ", vertex=" << vertex_id
				<< ", point=(" << V(vertex_id, 0)
				<< ", " << V(vertex_id, 1)
				<< ", " << V(vertex_id, 2) << ")"
				<< std::endl;
		}
		jud_error = true;
		current_remove_triangles.clear();
		current_slicer.clear();
		return;
	}

	// Rebuild the remaining face list in one pass. The selection mask is
	// indexed exactly like all_slicer.triangles, so no triangle search is needed.
	std::vector<TRiangle> remaining_triangles;
	remaining_triangles.reserve(all_slicer.triangles.size() - remove_triangles.size());
	for (int face_id = 0; face_id < static_cast<int>(all_slicer.triangles.size()); ++face_id) {
		if (!jud_triangle_have_been_added[face_id]) {
			remaining_triangles.push_back(all_slicer.triangles[face_id]);
		}
	}
	all_slicer.triangles.swap(remaining_triangles);

	//add cutting plane triangles
	/*all_slicer.triangles.insert(all_slicer.triangles.begin(),cutting_plane_points)
	cutting_plane_points*/

	cout << "cccc" << endl;
	//sort cutting_plane_points by adjacency relation
#if 0 // Replaced by ExtractAllClosedCutBoundaryLoops and AppendCutBoundaryCaps.
	vector<vector<int>> real_cutting_plane_triangles;	//为每个切割层存储排序后的切割平面顶点索引
	real_cutting_plane_triangles.resize(all_cut_layers.size());	//real_cutting_plane_triangles大小与all_cut_layers相同
	for (int i = 0; i < all_cut_layers.size(); i++) {
		if (all_cut_layers[i].empty() || cutting_plane_edges[i].empty()) {
			std::cout << "[HybridManufacturing::CutMesh] No cutting boundary edges for layer "
				<< i
				<< ", dependency=" << all_cut_layers_dependency_layer[i]
				<< std::endl;
			jud_error = true;
			continue;
		}

		cv::Point2d current_triangle_point;
		cv::Point2d current_layer_point;
		current_layer_point.x = all_cut_layers[i][0].x;
		current_layer_point.y = all_cut_layers[i][0].y;
		double min_dis = MAX_D;
		int index_start_point_id = -1;
		for (int j = 0; j < cutting_plane_edges[i].size(); j++) {
			current_triangle_point.x = all_slicer.positions[cutting_plane_edges[i][j].second][0];
			current_triangle_point.y = all_slicer.positions[cutting_plane_edges[i][j].second][1];
			if (distance2d(current_layer_point, current_triangle_point) < min_dis)
			{
				index_start_point_id = j;
				min_dis = distance2d(current_layer_point, current_triangle_point);
			}
		}
		if (index_start_point_id < 0) {
			std::cout << "[HybridManufacturing::CutMesh] Failed to select cutting boundary start edge for layer "
				<< i << std::endl;
			jud_error = true;
			continue;
		}
		int index_start_point = cutting_plane_edges[i][index_start_point_id].second;
		int index_of_last_edge = index_start_point_id;
		real_cutting_plane_triangles[i].push_back(cutting_plane_edges[i][index_start_point_id].first);
		//cout << "kkkk" << cutting_plane_edges[i].size() << endl;
		int jud_select_id = -1;
		int cont_segment = 0;
		while (index_start_point != cutting_plane_edges[i][index_start_point_id].first) {
			bool jud_select = false;
			cont_segment++;
			for (int j = 0; j < cutting_plane_edges[i].size(); j++) {
				if (j != index_of_last_edge && cutting_plane_edges[i][j].first == index_start_point) {
					real_cutting_plane_triangles[i].push_back(cutting_plane_edges[i][j].first);
					index_start_point = cutting_plane_edges[i][j].second;
					index_of_last_edge = j;
					jud_select = true;
					jud_select_id = 0;
					break;
				}
				else if (j != index_of_last_edge && cutting_plane_edges[i][j].second == index_start_point) {
					real_cutting_plane_triangles[i].push_back(cutting_plane_edges[i][j].second);
					index_start_point = cutting_plane_edges[i][j].first;
					index_of_last_edge = j;
					jud_select = true;
					jud_select_id = 1;
					break;
				}
			}
			if (cont_segment > 50000) {
				jud_select = false;
				break;
				//real_cutting_plane_triangles[i].clear();
			}
			if (cont_segment > 1000000) {
				jud_select = false;
				//real_cutting_plane_triangles[i].clear();
			}
			if (jud_select == false) {
				//cout << "no:" << i << "  ";
				jud_error = true;
				/*for (int j = 0; j < cutting_plane_edges[i].size(); j++) {
					if (j != index_of_last_edge) {
						cv::Point2d current_point_2, current_point_3, current_point_4;
						current_point_3.x = all_slicer.positions[cutting_plane_edges[i][j].first][0];
						current_point_3.y = all_slicer.positions[cutting_plane_edges[i][j].first][1];
						current_point_4.x = all_slicer.positions[cutting_plane_edges[i][j].second][0];
						current_point_4.y = all_slicer.positions[cutting_plane_edges[i][j].second][1];
						if (jud_select_id == 0) {
							current_point_2.x = all_slicer.positions[cutting_plane_edges[i][index_of_last_edge].second][0];
							current_point_2.y = all_slicer.positions[cutting_plane_edges[i][index_of_last_edge].second][1];
						}
						else if (jud_select_id == 1) {
							current_point_2.x = all_slicer.positions[cutting_plane_edges[i][index_of_last_edge].first][0];
							current_point_2.y = all_slicer.positions[cutting_plane_edges[i][index_of_last_edge].first][1];
						}
						if (distance_2d(current_point_2, current_point_3) < 0.00001) {
							real_cutting_plane_triangles[i].push_back(cutting_plane_edges[i][j].first);
							index_start_point = cutting_plane_edges[i][j].second;
							index_of_last_edge = j;
							jud_select = true;
							jud_select_id = 0;
							break;
						}
						else if (distance_2d(current_point_2, current_point_4) < 0.00001) {
							real_cutting_plane_triangles[i].push_back(cutting_plane_edges[i][j].second);
							index_start_point = cutting_plane_edges[i][j].first;
							index_of_last_edge = j;
							jud_select = true;
							jud_select_id = 1;
							break;
						}
					}
				}*/
				break;
			}
		}
		if (index_start_point != cutting_plane_edges[i][index_start_point_id].first) {
			std::cout << "[HybridManufacturing::CutMesh] Cutting boundary is not a closed loop for layer "
				<< i << std::endl;
			real_cutting_plane_triangles[i].clear();
			jud_error = true;
		}
	}
	/*if (jud_error[id_node] == true)
		return;*/
		/*cout << "mm" << endl;*/
		/*if (height_of_beam_search == 6)
			return;*/
	///////////////////////////terminate//////////////////////////////
	if (all_slicer.triangles.size() < terminate_threshold_of_number_of_faces) {
		jud_outer_beam_search_terminate = true;
	}
	//////////////////////////////////////////////////////////////////

	vector<int> all_cap_face_ids;
	vector<int> id_contact_faces;
	//if (height_of_beam_search != 2 || id_continue != 1)
	if (using_solid_model == true) {
		Anticlockwise(real_cutting_plane_triangles, all_slicer);
		std::map<vasco::contact_triangulation::QuantizedPointKey, int> point_index_map;
		for (const auto& tri : all_slicer.triangles) {
			for (int k = 0; k < 3; ++k) {
				point_index_map.emplace(
					vasco::contact_triangulation::MakeQuantizedKey(Point_3(all_slicer.positions[tri[k]][0], all_slicer.positions[tri[k]][1], all_slicer.positions[tri[k]][2])),
					tri[k]);
			}
		}
		for (int i = 0; i < real_cutting_plane_triangles.size(); i++) {
			if (real_cutting_plane_triangles[i].size() < 3)
				continue;
			if (i >= static_cast<int>(flag_cut_layers_is_hole.size()) || flag_cut_layers_is_hole[i] != -1)
				continue;
			const bool is_dependency_contact_face =
				i < static_cast<int>(all_cut_layers_dependency_layer.size())
				&& all_cut_layers_dependency_layer[i] > 0;
			std::vector<int> hole_ring;
			std::vector<std::vector<int>> hole_rings;
			for (int m = 0; m < real_cutting_plane_triangles.size(); m++) {
				if (m < static_cast<int>(flag_cut_layers_is_hole.size())
					&& flag_cut_layers_is_hole[m] == follow_index[i]
					&& real_cutting_plane_triangles[m].size() >= 3) {
					hole_ring = real_cutting_plane_triangles[m];
					hole_rings.push_back(hole_ring);
					break;
				}
			}

			//contact_face_triangulation_mode = ContactFaceTriangulationMode::Earcut;
			contact_face_triangulation_mode = ContactFaceTriangulationMode::CGALRemesh;
			if (contact_face_triangulation_mode == ContactFaceTriangulationMode::CGALRemesh) {
				std::vector<int> generated_cap_face_ids;
				vasco::contact_triangulation::TriangulateContactFacesWithCDT(
					real_cutting_plane_triangles[i],
					hole_rings,
					all_slicer,
					generated_cap_face_ids,
					point_index_map);
				all_cap_face_ids.insert(
					all_cap_face_ids.end(),
					generated_cap_face_ids.begin(),
					generated_cap_face_ids.end());
				if (is_dependency_contact_face) {
					id_contact_faces.insert(
						id_contact_faces.end(),
						generated_cap_face_ids.begin(),
						generated_cap_face_ids.end());
				}
				continue;
			}

			using Coord = double;
			using NN = uint32_t;
			using Point = std::array<Coord, 2>;
			std::vector<std::vector<Point>> polygon;
			polygon.clear();
			std::vector<Point> temp_vec;
			polygon.push_back(temp_vec);
			map<int, int> map_index_faces;
			map_index_faces.clear();
			for (int j = 0; j < real_cutting_plane_triangles[i].size(); j++) {
				polygon[0].push_back({ all_slicer.positions[real_cutting_plane_triangles[i][j]][0], all_slicer.positions[real_cutting_plane_triangles[i][j]][1] });
				map_index_faces.insert({ j,real_cutting_plane_triangles[i][j] });
			}
			for (int m = 0; m < real_cutting_plane_triangles.size(); m++) {
				if (m < static_cast<int>(flag_cut_layers_is_hole.size())
					&& flag_cut_layers_is_hole[m] == follow_index[i]
					&& real_cutting_plane_triangles[m].size() >= 3) {
					for (int j = 0; j < real_cutting_plane_triangles[m].size(); j++) {
						polygon.push_back(temp_vec);
						polygon[1].push_back({ all_slicer.positions[real_cutting_plane_triangles[m][j]][0], all_slicer.positions[real_cutting_plane_triangles[m][j]][1] });
						map_index_faces.insert({ polygon[0].size() + j,real_cutting_plane_triangles[m][j] });
					}
					break;
				}
			}

			std::vector<std::vector<Eigen::Vector3d>> polygons_v3d;
			for (int j = 0; j < polygon.size(); j++) {
				std::vector<Eigen::Vector3d> temp_polygon_v3d;
				for (int k = 0; k < polygon[j].size(); k++) {
					Eigen::Vector3d temp_point_v3d(polygon[j][k][0], polygon[j][k][1], all_slicer.positions[real_cutting_plane_triangles[i][0]][2]);
					temp_polygon_v3d.push_back(temp_point_v3d);
				}
				polygons_v3d.push_back(temp_polygon_v3d);
			}
			vasco::io::outputContourToObj(polygons_v3d, Eigen::Vector3d(0.5, 0.5, 0.0), file_name + "-contact_face_polygon-" + std::to_string(height_of_beam_search) + "_" + std::to_string(i) + ".obj");

			std::vector<NN> indices = mapbox::earcut<NN>(polygon);
			for (int j = 0; j < indices.size();) {
				TRiangle the_new_cutting_plane_triangle;
				the_new_cutting_plane_triangle[0] = map_index_faces[indices[j]]; j++;
				the_new_cutting_plane_triangle[1] = map_index_faces[indices[j]]; j++;
				the_new_cutting_plane_triangle[2] = map_index_faces[indices[j]]; j++;
				all_slicer.triangles.insert(all_slicer.triangles.end(), the_new_cutting_plane_triangle);
				const int cap_face_id = static_cast<int>(all_slicer.triangles.size() - 1);
				all_cap_face_ids.push_back(cap_face_id);
				if (is_dependency_contact_face) {
					id_contact_faces.push_back(cap_face_id);
				}
			}

		}
		//}
	}

#endif
	if (all_slicer.triangles.size() < terminate_threshold_of_number_of_faces) {
		jud_outer_beam_search_terminate = true;
		std::cout << "[HybridManufacturing::CutMesh] Terminate beam search by remaining face count: "
			<< all_slicer.triangles.size()
			<< " < " << terminate_threshold_of_number_of_faces
			<< std::endl;
	}

	bool all_cut_loops_valid = true;
	const std::vector<std::vector<std::vector<int>>> extracted_cut_boundary_loops =
		ExtractAllClosedCutBoundaryLoops(cutting_plane_edges, all_cut_loops_valid);
	RemainingCutBoundaryLoopStats remaining_loop_stats;
	const std::vector<std::vector<std::vector<int>>> cutting_boundary_loops =
		FilterCutLoopsToRemainingMeshBoundary(
			all_slicer,
			extracted_cut_boundary_loops,
			remaining_loop_stats,
			all_cut_loops_valid);
	if (!all_cut_loops_valid) {
		jud_error = true;
	}

	std::size_t cut_loop_count = 0;
	for (const auto& layer_loops : cutting_boundary_loops) {
		cut_loop_count += layer_loops.size();
	}
	std::cout << "[HybridManufacturing::CutMesh] Remaining cut boundary loop filter: "
		<< "extracted=" << remaining_loop_stats.input_loop_count
		<< ", retained=" << cut_loop_count
		<< ", discarded_removed_side=" << remaining_loop_stats.discarded_orphan_loop_count
		<< ", partial=" << remaining_loop_stats.partial_boundary_loop_count
		<< std::endl;

	vector<int> all_cap_face_ids;
	vector<int> id_contact_faces;
	if (using_solid_model) {
		contact_face_triangulation_mode = ContactFaceTriangulationMode::CGALRemesh;
		AppendCutBoundaryCaps(
			cutting_boundary_loops,
			flag_cut_layers_is_hole,
			follow_index,
			all_cut_layers_dependency_layer,
			all_slicer,
			all_cap_face_ids,
			id_contact_faces);
	}




	/*if (height_of_beam_search == 1) {
		for (int i = 0; i < real_cutting_plane_triangles.size(); i++) {
			VEctor new_point;
			double temp_x = 0, temp_y = 0;
			for (int j = 0; j < real_cutting_plane_triangles[i].size(); j++) {
				temp_x += all_slicer.positions[real_cutting_plane_triangles[i][j]][0];
				temp_y += all_slicer.positions[real_cutting_plane_triangles[i][j]][1];
			}
			new_point[0] = temp_x / real_cutting_plane_triangles[i].size();
			new_point[1] = temp_y / real_cutting_plane_triangles[i].size();
			new_point[2] = all_slicer.positions[real_cutting_plane_triangles[i][0]][2];
			all_slicer.positions.insert(all_slicer.positions.end(), new_point);
			for (int j = 0; j < real_cutting_plane_triangles[i].size(); j++) {
				TRiangle the_new_cutting_plane_triangle;
				the_new_cutting_plane_triangle[0] = all_slicer.positions.size() - 1;
				the_new_cutting_plane_triangle[1] = real_cutting_plane_triangles[i][j];
				the_new_cutting_plane_triangle[2] = real_cutting_plane_triangles[i][(j + 1) % real_cutting_plane_triangles[i].size()];
				all_slicer.triangles.insert(all_slicer.triangles.end(), the_new_cutting_plane_triangle);
			}
		}
	}*/

	std::cout << "[HybridManufacturing::CutMesh] Cap face classification: "
		<< "all_caps=" << all_cap_face_ids.size()
		<< ", dependency_contact_faces=" << id_contact_faces.size()
		<< ", zero_dependency_caps=" << (all_cap_face_ids.size() - id_contact_faces.size())
		<< std::endl;

	SaveRemovedSolidAfterCut(
		all_slicer,
		remove_triangles,
		all_cap_face_ids,
		vector_after,
		vectorBefore,
		height_of_beam_search,
		cont_number_of_queue,
		judge_continue_additive,
		id_continue);

	Slicer_2 all_slicer_2 = all_slicer;
	all_slicer_2.triangles = remove_triangles;

	RotateSlicerPositions(all_slicer_2, vector_after, vectorBefore);


	if (judge_continue_additive == false)
		all_slicer_2.save(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_current" + ".obj");
	else
		all_slicer_2.save(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_" + to_string(id_continue) + "_current" + "_subblock.obj");

	//std::vector<double> cut_plane_z_values;
	//cut_plane_z_values.reserve(all_cut_layers.size());
	//for (const auto& cut_layer : all_cut_layers) {
	//	if (!cut_layer.empty()) {
	//		cut_plane_z_values.push_back(cut_layer[0].z);
	//	}
	//}
	//vasco::slicer_mesh_adapter::RemeshCutPatchBeforeSaving(
	//	all_slicer,
	//	id_contact_faces,
	//	cut_plane_z_values,
	//	std::max(dh, 5.0));


	all_slicer.save(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_before_rotate_back.obj");

	RotateSlicerPositions(all_slicer, vector_after, vectorBefore);

	constexpr double remaining_mesh_z_extent_threshold = 3.0;
	ReferencedSlicerZBounds referenced_z_bounds;
	const bool has_referenced_mesh_bounds = ComputeReferencedSlicerZBounds(all_slicer, referenced_z_bounds);
	if (!has_referenced_mesh_bounds) {
		jud_outer_beam_search_terminate = true;
		std::cout << "[HybridManufacturing::CutMesh] Terminate beam search because the remaining mesh "
			<< "has no finite referenced vertices: faces=" << all_slicer.triangles.size()
			<< ", invalid_vertex_indices=" << referenced_z_bounds.invalid_vertex_index_count
			<< std::endl;
	}
	else {
		const double referenced_z_extent = referenced_z_bounds.max_z - referenced_z_bounds.min_z;
		std::cout << "[HybridManufacturing::CutMesh] Remaining referenced mesh z bounds: min="
			<< referenced_z_bounds.min_z
			<< ", max=" << referenced_z_bounds.max_z
			<< ", extent=" << referenced_z_extent
			<< ", threshold=" << remaining_mesh_z_extent_threshold
			<< ", referenced_vertices=" << referenced_z_bounds.referenced_vertex_count
			<< ", stored_vertices=" << all_slicer.positions.size()
			<< ", invalid_vertex_indices=" << referenced_z_bounds.invalid_vertex_index_count
			<< std::endl;

		if (referenced_z_extent < remaining_mesh_z_extent_threshold) {
			jud_outer_beam_search_terminate = true;
			std::cout << "[HybridManufacturing::CutMesh] Terminate beam search by remaining referenced mesh z extent."
				<< std::endl;
		}
	}


	Eigen::MatrixXd temp_V(all_slicer.positions.size(), 3);
	Eigen::MatrixXi temp_F(all_slicer.triangles.size(), 3);
	for (int i = 0; i < all_slicer.positions.size(); i++) {
		temp_V.row(i).x() = all_slicer.positions[i][0];
		temp_V.row(i).y() = all_slicer.positions[i][1];
		temp_V.row(i).z() = all_slicer.positions[i][2];
	}
	for (int i = 0; i < all_slicer.triangles.size(); i++) {
		temp_F.row(i).x() = all_slicer.triangles[i][0];
		temp_F.row(i).y() = all_slicer.triangles[i][1];
		temp_F.row(i).z() = all_slicer.triangles[i][2];
	}
	if (judge_continue_additive == false) {
		Katana::Instance().stl.saveStlFromObj(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + ".stl", temp_V, temp_F);
		vector<std::array<double, 3>> V3;
		vector<std::array<int, 3>> F3;
		vector<std::array<double, 3>> N3;
		ifstream ifs(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + ".stl");
		igl::read_stl_ascii(ifs, V3, F3, N3);
		Eigen::MatrixXd V4(V3.size(), 3);
		Eigen::MatrixXi F4(F3.size(), 3);
		Eigen::MatrixXd N4(N3.size(), 3);
		for (int i = 0; i < V3.size(); i++)
			for (int j = 0; j < 3; j++)
				V4(i, j) = V3[i][j];
		for (int i = 0; i < F3.size(); i++)
			for (int j = 0; j < 3; j++)
				F4(i, j) = F3[i][j];
		igl::writeSTL(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_B.stl", V4, F4, igl::FileEncoding::Binary);

		const std::string output_obj =
			file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + ".obj";
		if (!igl::writeOBJ(output_obj, temp_V, temp_F)) {
			std::cerr << "[HybridManufacturing::CutMesh] Failed to write full-precision OBJ: "
				<< output_obj << std::endl;
			jud_error = true;
		}
	}
	else {
		Katana::Instance().stl.saveStlFromObj(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_" + to_string(id_continue) + "_subblock.stl", temp_V, temp_F);
		vector<std::array<double, 3>> V3;
		vector<std::array<int, 3>> F3;
		vector<std::array<double, 3>> N3;
		ifstream ifs(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_" + to_string(id_continue) + "_subblock.stl");
		igl::read_stl_ascii(ifs, V3, F3, N3);
		Eigen::MatrixXd V4(V3.size(), 3);
		Eigen::MatrixXi F4(F3.size(), 3);
		Eigen::MatrixXd N4(N3.size(), 3);
		for (int i = 0; i < V3.size(); i++)
			for (int j = 0; j < 3; j++)
				V4(i, j) = V3[i][j];
		for (int i = 0; i < F3.size(); i++)
			for (int j = 0; j < 3; j++)
				F4(i, j) = F3[i][j];
		igl::writeSTL(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_" + to_string(id_continue) + "_subblock_B.stl", V4, F4, igl::FileEncoding::Binary);

		const std::string output_obj =
			file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue)
			+ "_" + to_string(id_continue) + "_subblock.obj";
		if (!igl::writeOBJ(output_obj, temp_V, temp_F)) {
			std::cerr << "[HybridManufacturing::CutMesh] Failed to write full-precision OBJ: "
				<< output_obj << std::endl;
			jud_error = true;
		}
	}

	string str_contact_faces;
	if (judge_continue_additive == false)
		str_contact_faces = file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_contact_faces.txt";
	else
		str_contact_faces = file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_" + to_string(id_continue) + "_subblock_contact_faces.txt";

	Visual vis_2;
	cv::Point3d input_ori(vector_after.x(), vector_after.y(), vector_after.z());
	if (judge_continue_additive == false)
		vis_2.generateArrows(input_ori, file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue));
	else
		vis_2.generateArrows(input_ori, file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_" + to_string(id_continue) + "_subblock");

	ofstream ofile_contact_faces(str_contact_faces);
	ofile_contact_faces << id_contact_faces.size() << endl;
	for (int i = 0; i < id_contact_faces.size(); i++)
		ofile_contact_faces << id_contact_faces[i] << endl;
	cout << "d" << endl;
	/////////////////////////
	return;
}

int HybridManufacturing::subtractive_accessibility_decomposition(
	vector<TRiangle> need_detect_triangle,
	int height_of_beam_search,
	int cont_number_of_queue,
	cutter cutting_tool,
	Slicer_2 current_slicer)
{
	cutting_tool.cylinder_r = 1.5;
	cutting_tool.cylinder_height = 28;
	cutting_tool.ball_r = 1.5;
	cutting_tool.carriage_r = 23;
	cutting_tool.carriage_height = 33;

	/////show all accessible points in every orientation/////
	Eigen::MatrixXd V_B2;
	Eigen::MatrixXi F_B2;
	igl::readOBJ("ball.obj", V_B2, F_B2);
	const Slicer_2 slicer = current_slicer; //加个const看看有没有涉及更改
	/*slicer.load((file_name + "-" + to_string(height_of_beam_search - 1) + "_" + to_string(cont_number_of_queue) + ".obj").c_str());*/

	/////////////////Get sample points//////////////////
	if (need_detect_triangle.empty()) {
		std::cout << "[HybridManufacturing::subtractive_accessibility_decomposition] no triangles to decompose: "
			<< "height=" << height_of_beam_search
			<< ", queue=" << cont_number_of_queue
			<< std::endl;
		return 0;
	}

	vector<cv::Point3d> all_sample_points_in_triangles(need_detect_triangle.size());
	for (int i = 0; i < need_detect_triangle.size(); i++) {

		cv::Point3d Point1(slicer.positions[need_detect_triangle[i][0]][0], slicer.positions[need_detect_triangle[i][0]][1], slicer.positions[need_detect_triangle[i][0]][2]);
		cv::Point3d Point2(slicer.positions[need_detect_triangle[i][1]][0], slicer.positions[need_detect_triangle[i][1]][1], slicer.positions[need_detect_triangle[i][1]][2]);
		cv::Point3d Point3(slicer.positions[need_detect_triangle[i][2]][0], slicer.positions[need_detect_triangle[i][2]][1], slicer.positions[need_detect_triangle[i][2]][2]);
		double a = distance3d(Point1, Point2);
		double b = distance3d(Point1, Point3);
		double c = distance3d(Point2, Point3);
		cv::Point3d V_incentre(((a * Point1.x + b * Point2.x + c * Point3.x) / (a + b + c)), ((a * Point1.y + b * Point2.y + c * Point3.y) / (a + b + c)), ((a * Point1.z + b * Point2.z + c * Point3.z) / (a + b + c)));
		all_sample_points_in_triangles[i] = V_incentre;
	}
	///////////////////////////////////////////////////
	cout << "% " << slicer.triangles.size() << endl;
	cout << "% " << need_detect_triangle.size() << endl;
	vector<vector<int>> accessible_ori_of_need_detect_V(need_detect_triangle.size(), vector<int>(sampling_subtractive.sample_points.size(), 10000000));
	vector<vector<Eigen::MatrixXd>> vis_points(1);
	vector<vector<vector<Eigen::Vector3d>>> vis_lines(1);


	cout << "%sampling_subtractive.sample_points.size() " << sampling_subtractive.sample_points.size() << endl;
	for (int ori = 0; ori < sampling_subtractive.sample_points.size(); ori++) {
		//先去除水平面以下的方向  //注意只有最底下的block需要限制!!!!!!!!!!!!!!
		/*if (sampling_subtractive.sample_points[ori].z < 0.2)
			continue;*/
		vector<std::vector<Eigen::Vector3d>> temp_new_V_remain(slicer.triangles.size());
		for (int i = 0; i < slicer.triangles.size(); i++) {
			temp_new_V_remain[i].resize(3);
			for (int k = 0; k < 3; k++) {
				temp_new_V_remain[i][k](0, 0) = slicer.positions[slicer.triangles[i][k]][0];
				temp_new_V_remain[i][k](1, 0) = slicer.positions[slicer.triangles[i][k]][1];
				temp_new_V_remain[i][k](2, 0) = slicer.positions[slicer.triangles[i][k]][2];
			}
		}

		std::vector<Eigen::Vector3d> temp_V_need_detect(need_detect_triangle.size());
		for (int i = 0; i < need_detect_triangle.size(); i++) {
			temp_V_need_detect[i](0, 0) = all_sample_points_in_triangles[i].x;
			temp_V_need_detect[i](1, 0) = all_sample_points_in_triangles[i].y;
			temp_V_need_detect[i](2, 0) = all_sample_points_in_triangles[i].z;
		}

		vector<std::vector<Eigen::Vector3d>> temp_V_need_detect_triangle(need_detect_triangle.size());
		for (int i = 0; i < need_detect_triangle.size(); i++) {
			temp_V_need_detect_triangle[i].resize(3);
			for (int k = 0; k < 3; k++) {
				temp_V_need_detect_triangle[i][k](0, 0) = slicer.positions[need_detect_triangle[i][k]][0];
				temp_V_need_detect_triangle[i][k](1, 0) = slicer.positions[need_detect_triangle[i][k]][1];
				temp_V_need_detect_triangle[i][k](2, 0) = slicer.positions[need_detect_triangle[i][k]][2];
			}
		}

		Eigen::Matrix3d rotMatrix;
		Eigen::Vector3d vectorBefore(0, 0, 1);
		Eigen::Vector3d vectorAfter(sampling_subtractive.sample_points[ori]);
		vectorAfter.normalize();
		rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vectorAfter).toRotationMatrix();
		for (int i = 0; i < temp_new_V_remain.size(); i++)
			for (int j = 0; j < temp_new_V_remain[i].size(); j++)
				temp_new_V_remain[i][j] = rotMatrix.inverse() * temp_new_V_remain[i][j];
		for (int i = 0; i < temp_V_need_detect.size(); i++)
			temp_V_need_detect[i] = rotMatrix.inverse() * temp_V_need_detect[i];
		for (int i = 0; i < temp_V_need_detect_triangle.size(); i++)
			for (int j = 0; j < temp_V_need_detect_triangle[i].size(); j++)
				temp_V_need_detect_triangle[i][j] = rotMatrix.inverse() * temp_V_need_detect_triangle[i][j];

		vector<Eigen::Vector3d> all_normals_of_need_detect_triangle;
		all_normals_of_need_detect_triangle.clear();
		for (int i = 0; i < need_detect_triangle.size(); i++) {
			cv::Point3d Point1(temp_V_need_detect_triangle[i][0](0, 0), temp_V_need_detect_triangle[i][0](1, 0), temp_V_need_detect_triangle[i][0](2, 0));
			cv::Point3d Point2(temp_V_need_detect_triangle[i][1](0, 0), temp_V_need_detect_triangle[i][1](1, 0), temp_V_need_detect_triangle[i][1](2, 0));
			cv::Point3d Point3(temp_V_need_detect_triangle[i][2](0, 0), temp_V_need_detect_triangle[i][2](1, 0), temp_V_need_detect_triangle[i][2](2, 0));
			double na = (Point2.y - Point1.y) * (Point3.z - Point1.z) - (Point2.z - Point1.z) * (Point3.y - Point1.y);
			double nb = (Point2.z - Point1.z) * (Point3.x - Point1.x) - (Point2.x - Point1.x) * (Point3.z - Point1.z);
			double nc = (Point2.x - Point1.x) * (Point3.y - Point1.y) - (Point2.y - Point1.y) * (Point3.x - Point1.x);
			Eigen::Vector3d vn(na, nb, nc);
			vn.normalize();
			all_normals_of_need_detect_triangle.push_back(vn);
			//cout << vn.x() << " "<<vn.y()<<" "<<vn.z() << endl;
		}

		vector<double> max_z_of_triangles(temp_new_V_remain.size());
		for (int i = 0; i < temp_new_V_remain.size(); i++) {
			max_z_of_triangles[i] = MIN_D;
			for (int j = 0; j < 3; j++)
				max_z_of_triangles[i] = max(max_z_of_triangles[i], temp_new_V_remain[i][j](2, 0));
		}

		///////////////////collision detection////////////////////////
		PrepareToolForCollision(cutting_tool);

		for (int i = 0; i < temp_V_need_detect.size(); i++) {
			bool jud_collision = false;
			Eigen::Vector3d center_point;

			center_point.x() = temp_V_need_detect[i](0, 0) + (cutting_tool.cylinder_r) * all_normals_of_need_detect_triangle[i].x();
			center_point.y() = temp_V_need_detect[i](1, 0) + (cutting_tool.cylinder_r) * all_normals_of_need_detect_triangle[i].y();
			center_point.z() = temp_V_need_detect[i](2, 0) + (cutting_tool.cylinder_r) * all_normals_of_need_detect_triangle[i].z();

			for (int ii = 0; ii < temp_new_V_remain.size(); ii++) {

				if (CheckToolCollisionWithCell(center_point, temp_new_V_remain[ii], max_z_of_triangles[ii], cutting_tool, 30.0, 3.0)) {
					jud_collision = true;
					break;
				}
			}

			if (jud_collision == false) {
				accessible_ori_of_need_detect_V[i][ori] = 0;
			}
		}


		//visualize
		/*rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorAfter,vectorBefore).toRotationMatrix();
		for (int i = 0; i < vis_points.size(); i++)
			vis_points[i] = rotMatrix.inverse() * vis_points[i];
		ofstream all_balls(".\\vis\\coral_accessible_points_in_ori-" + to_string(ori) + ".obj");
		for (int i = 0; i < vis_points.size(); i++) {
			for (int j = 0; j < V_B2.rows(); j++)
				all_balls << "v " << V_B2(j, 0) + vis_points[i](0, 0) << " " << V_B2(j, 1) + vis_points[i](1, 0) << " " << V_B2(j, 2) + vis_points[i](2, 0) << " 0.9" << " 0.05" << " 0.05" << endl;
			for (int j = 0; j < F_B2.rows(); j++)
				all_balls << "f " << F_B2(j, 0) + i * V_B2.rows() + 1 << " " << F_B2(j, 1) + i * V_B2.rows() + 1 << " " << F_B2(j, 2) + i * V_B2.rows() + 1 << endl;
		}
		all_balls.close();
		Visual vis;
		vis.generateModelForRendering_7(sampling_subtractive.sample_points[ori], ".\\vis\\coral_accessible_points_in_ori-" + to_string(ori) + "_orientation.obj");*/
	}


	//////////////////////////////////graph cut////////////////////////////////////////
	//不应该出现cont_ori == 0的情况
	int cont_revise = 0;
	for (int i = 0; i < accessible_ori_of_need_detect_V.size(); i++) {
		int cont_ori = 0;
		for (int j = 0; j < accessible_ori_of_need_detect_V[i].size(); j++) {
			if (accessible_ori_of_need_detect_V[i][j] == 0)
				cont_ori++;
		}
		//暂时先强制任意方向可达
		if (cont_ori == 0) {
			cont_revise++;
			//cout << "*** " << need_detect_triangle[i][0] << endl;
			for (int j = 0; j < accessible_ori_of_need_detect_V[i].size(); j++)
				accessible_ori_of_need_detect_V[i][j] = 0;
		}
		/*if (accessible_ori_of_need_detect_V[i][123] != 0)
			cout << "no" << endl;*/
	}
	cout << cont_revise << endl;


	vector<std::vector<int>> pixels_relations;
	vector<int> length_edges;
	for (int i = 0; i < need_detect_triangle.size(); i++)
		for (int j = 0; j < need_detect_triangle.size(); j++) {
			bool jud_adjacent = false;
			vector<int> temp_edge(2);
			for (int ii = 0; ii < 3; ii++) {
				for (int jj = 0; jj < 3; jj++) {
					if (need_detect_triangle[i][ii] == need_detect_triangle[j][jj]) {
						/*temp_edge[0] = i;
						temp_edge[1] = j;
						pixels_relations.push_back(temp_edge);
						break;*/
						if (need_detect_triangle[i][(ii + 1) % 3] == need_detect_triangle[j][(jj + 1) % 3] || need_detect_triangle[i][(ii + 1) % 3] == need_detect_triangle[j][(jj + 2) % 3]) {
							temp_edge[0] = i;
							temp_edge[1] = j;
							pixels_relations.push_back(temp_edge);
							float temp_length = distanceVec3(slicer.positions[need_detect_triangle[i][ii]], slicer.positions[need_detect_triangle[i][(ii + 1) % 3]]) * 100;
							length_edges.push_back(temp_length);
							jud_adjacent = true;
							break;
						}
						else if (need_detect_triangle[i][(ii + 2) % 3] == need_detect_triangle[j][(jj + 1) % 3] || need_detect_triangle[i][(ii + 2) % 3] == need_detect_triangle[j][(jj + 2) % 3]) {
							temp_edge[0] = i;
							temp_edge[1] = j;
							pixels_relations.push_back(temp_edge);
							float temp_length = distanceVec3(slicer.positions[need_detect_triangle[i][ii]], slicer.positions[need_detect_triangle[i][(ii + 2) % 3]]) * 100;
							length_edges.push_back(temp_length);
							jud_adjacent = true;
							break;
						}
						else {
							jud_adjacent = true;
							break;
						}
					}
				}
				if (jud_adjacent == true)
					break;
			}
			/*if (jud_adjacent == false) {
				temp_edge[0] = i;
				temp_edge[1] = j;
				pixels_relations.push_back(temp_edge);
				length_edges.push_back(100000);
			}*/
		}


	cout << endl << "Graph Cut................" << endl;
	cout << "% " << need_detect_triangle.size() << endl;
	cout << "% " << sampling_subtractive.sample_points.size() << endl;
	vector<int> result = GeneralGraph_DArraySArraySpatVarying(need_detect_triangle.size(), sampling_subtractive.sample_points.size(), accessible_ori_of_need_detect_V, pixels_relations, length_edges);

	for (int i = 0; i < need_detect_triangle.size(); i++)
		for (int j = i + 1; j < need_detect_triangle.size(); j++)
			if (result[i] > result[j]) {
				swap(result[i], result[j]);
				swap(need_detect_triangle[i], need_detect_triangle[j]);
			}

	//cout << "normal:" << result[0]<<endl;
	ofstream ofile(PatchVisPath("normal_of_pathches.txt"));
	vector<vector<vasco::core::Tri3>> vis_triangles(1);
	vector<vasco::core::Vec3> vis_positions = slicer.positions;

	int cont_patch = 0;
	vector<Eigen::Vector3d> points_in_cell;
	vector<Eigen::Vector3d> normals;
	points_in_cell.push_back(Eigen::Vector3d(slicer.positions[need_detect_triangle[0][0]][0], slicer.positions[need_detect_triangle[0][0]][1], slicer.positions[need_detect_triangle[0][0]][2])); //point_in_cell初始放入第一个点
	normals.push_back(sampling_subtractive.sample_points[result[0]]);
	ofile << normals[0].x() << " " << normals[0].y() << " " << normals[0].z() << endl;


	for (int i = 0; i < need_detect_triangle.size(); i++) {
		if (i != 0 && result[i - 1] < result[i]) {
			cont_patch++;
			vector<Eigen::MatrixXd> temp_vecc(0);
			vis_points.push_back(temp_vecc);
			vector<vector<Eigen::Vector3d>> temp_vec(0);
			vis_lines.push_back(temp_vec);
			vector<vasco::core::Tri3> temp_tri(0);
			vis_triangles.push_back(temp_tri);
			points_in_cell.push_back(Eigen::Vector3d(slicer.positions[need_detect_triangle[i][0]][0], slicer.positions[need_detect_triangle[i][0]][1], slicer.positions[need_detect_triangle[i][0]][2]));
			normals.push_back(sampling_subtractive.sample_points[result[i]]);
			ofile << normals[normals.size() - 1].x() << " " << normals[normals.size() - 1].y() << " " << normals[normals.size() - 1].z() << endl;
		}

		//Eigen::MatrixXd temp_point(3, 1);
		//temp_point(0, 0) = V(index_V_need_detect[i], 0);
		//temp_point(1, 0) = V(index_V_need_detect[i], 1);
		//temp_point(2, 0) = V(index_V_need_detect[i], 2);
		////cout << vis_points[cont_patch].size() << endl;
		//vis_points[cont_patch].push_back(temp_point);
		//cout << "t" << endl;
		vector<Eigen::Vector3d> temp_vec;

		temp_vec.clear();
		vasco::core::Tri3 temp_tri;
		//cout << index_V_need_detect[i] << endl;
		for (int j = 0; j < 3; j++) {
			//cout << "f" << endl;
			Eigen::Vector3d temp_vec_2(slicer.positions[need_detect_triangle[i][j]][0], slicer.positions[need_detect_triangle[i][j]][1], slicer.positions[need_detect_triangle[i][j]][2]);
			//cout << "g" << endl;
			temp_vec.push_back(temp_vec_2);
			temp_tri[j] = need_detect_triangle[i][j];
		}
		//cout << "d" << endl;
		//cout << vis_lines.size() << " " << cont_patch << endl;
		vis_lines[cont_patch].push_back(temp_vec);
		vis_triangles[cont_patch].push_back(temp_tri);
		//cout << vis_lines[cont_patch].size() << endl;
		//cout << "e" << endl;
	}

	//visualize	
	/*ofstream all_balls(".\\vis\\coral_subtractive_decompose.obj");
	int cont_v = 0;
	for (int i = 0; i < vis_points.size(); i++) {
		double r = rand() / double(RAND_MAX);
		double g = rand() / double(RAND_MAX);
		double b = rand() / double(RAND_MAX);
		for (int j = 0; j < vis_points[i].size(); j++) {
			for (int k = 0; k < V_B2.rows(); k++)
				all_balls << "v " << V_B2(k, 0) + vis_points[i][j](0, 0) << " " << V_B2(k, 1) + vis_points[i][j](1, 0) << " " << V_B2(k, 2) + vis_points[i][j](2, 0) << " " << r << " " << g << " " << b << endl;
			for (int k = 0; k < F_B2.rows(); k++)
				all_balls << "f " << F_B2(k, 0) + cont_v * V_B2.rows() + 1 << " " << F_B2(k, 1) + cont_v * V_B2.rows() + 1 << " " << F_B2(k, 2) + cont_v * V_B2.rows() + 1 << endl;
			cont_v++;
		}
	}
	all_balls.close();*/
	for (int t = 0; t < vis_lines.size(); t++) {
		std::ofstream dstream(PatchVisPath(
			"patch-" + to_string(height_of_beam_search)
			+ "_" + to_string(cont_number_of_queue)
			+ "_" + to_string(t)
			+ ".stl"));
		if (!dstream.is_open()) {
			std::cout << "can not open " << std::endl;
			return static_cast<int>(vis_lines.size());
		}
		dstream << "solid STL generated by MeshLab" << std::endl;
		for (int i = 0; i < vis_lines[t].size(); i++) {
			dstream << "  facet normal " << "0 0 0" << std::endl;
			dstream << "    outer loop" << std::endl;
			for (int j = 0; j < 3; j++) {
				dstream << "      vertex  " << vis_lines[t][i][j][0] << " " << vis_lines[t][i][j][1] << " " << vis_lines[t][i][j][2] << std::endl;
			}
			dstream << "    endloop" << std::endl;
			dstream << "  endfacet" << std::endl;
		}
		dstream << "endsolid vcg" << std::endl;
		dstream.close();
	}
	//for (int t = 0; t < vis_triangles.size(); t++) {
	//	std::string mesh_name = "subtractive_patch_" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_" + to_string(t);
	//	polyscope::registerSurfaceMesh(mesh_name, vis_positions, vis_triangles[t]);
	//}
	//polyscope::show();
	ofile.close();


	cout << cont_patch << endl;
	vector<vector<double>> colors(vis_lines.size(), vector<double>(3));
	Visual vis;
	vis.generateModelForRendering_10(
		vis_lines,
		colors,
		PatchVisPath("coral_subtractive_patch_voronoi_cell.obj"));
	Visual vis_2;
	vis_2.generateModelForRendering_11(
		points_in_cell,
		normals,
		colors,
		PatchVisPath("coral_subtractive_patch_normal.obj"));
	//Visual vis;
	//vis.generateModelForRendering_7(sampling_subtractive.sample_points[ori], ".\\vis\\coral_accessible_points_in_ori-" + to_string(ori) + "_orientation.obj");
	//////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////
	return static_cast<int>(vis_lines.size());
}

void HybridManufacturing::subtractive_accessibility_decomposition_global(int height_of_beam_search, cutter cutting_tool)
{
	RunSubtractiveAccessibilityDecompositionGlobal(
		height_of_beam_search,
		cutting_tool,
		false,
		0,
		0,
		1);
}

void HybridManufacturing::subtractive_accessibility_decomposition_local_initialized_global(
	int height_of_beam_search,
	cutter cutting_tool)
{
	RunSubtractiveAccessibilityDecompositionGlobal(
		height_of_beam_search,
		cutting_tool,
		true,
		0,
		0,
		1);
}

void HybridManufacturing::subtractive_accessibility_decomposition_random_multistart_global(
	int height_of_beam_search,
	cutter cutting_tool,
	int restart_count,
	unsigned int random_seed,
	int thread_count)
{
	RunSubtractiveAccessibilityDecompositionGlobal(
		height_of_beam_search,
		cutting_tool,
		false,
		std::max(1, restart_count),
		random_seed,
		thread_count);
}

void HybridManufacturing::RunSubtractiveAccessibilityDecompositionGlobal(
	int height_of_beam_search,
	cutter cutting_tool,
	bool initialize_with_local_labels,
	int random_restart_count,
	unsigned int random_seed,
	int random_thread_count)
{
	ResetPolyscopeSceneForGlobalDecomposition();
	const bool use_random_multistart = random_restart_count > 0;
	const std::string output_prefix = initialize_with_local_labels
		? "local_initialized_global_"
		: (use_random_multistart ? "random_multistart_global_" : "");

	sampling_subtractive.OrientationSamplePoints();	//sampling_subtractive生成球面采样点
	cutting_tool.cylinder_r = 1.5;
	cutting_tool.cylinder_height = 28;
	cutting_tool.ball_r = 1.5;
	cutting_tool.carriage_r = 23;
	cutting_tool.carriage_height = 33;

	/////show all accessible points in every orientation/////
	Eigen::MatrixXd V_B2;
	Eigen::MatrixXi F_B2;
	igl::readOBJ("ball.obj", V_B2, F_B2);

	Slicer_2 slicer_load_current_patch;

	Slicer_2 slicer_load_next_patch;

	const std::string ancestor_source_report_file =
		AncestorSourceVisPath(
			"beam_search_last_layer_ancestor_sources_"
			+ MakeAccessibilityDebugFileToken(file_name)
			+ ".txt");
	const std::vector<std::string> ancestor_patch_files =
		LoadFirstFinalNodeAncestorPatchFiles(
			ancestor_source_report_file,
			file_name,
			PatchVisDir());
	const int graph_cut_block_count = ancestor_patch_files.empty()
		? height_of_beam_search
		: static_cast<int>(ancestor_patch_files.size());

	if (!ancestor_patch_files.empty()) {
		std::cout << "[Info] global decomposition uses ancestor patches from first final node in "
			<< ancestor_source_report_file << std::endl;
		for (int i = 0; i < static_cast<int>(ancestor_patch_files.size()); ++i) {
			std::cout << "  ancestor_patch[" << (i + 1) << "]: "
				<< ancestor_patch_files[i] << std::endl;
		}
	}

	// 读取并拼接最后一层第0号节点及其祖先 patch；若报告不存在则退回旧的 block_patch-<i>_0.obj
	std::vector<int> merged_vertex_source_patch_id;
	std::vector<int> merged_face_source_patch_id;

	Slicer_2 slicer_load_merged_patch = ancestor_patch_files.empty()
		? MergeBlockPatchesWithDedup(
			height_of_beam_search,
			merged_vertex_source_patch_id,
			merged_face_source_patch_id,
			1e-3)
		: MergeBlockPatchesWithDedup(
			ancestor_patch_files,
			merged_vertex_source_patch_id,
			merged_face_source_patch_id,
			1e-3);
	if (!StitchMergedPatchBoundariesForGraphCut(
		slicer_load_merged_patch,
		merged_face_source_patch_id,
		"global")) {
		return;
	}
	slicer_load_merged_patch.save(
		PatchVisPath(output_prefix + "block_patch_merged_removedup.obj"));
	std::cout << "[Info] merged patch mesh: V=" << slicer_load_merged_patch.positions.size()
		<< ", F=" << slicer_load_merged_patch.triangles.size() << std::endl;

	// 计算 merged patch 上每个三角面在每个采样方向下的最大碰撞 patch_index，以此得到可行patch_index的区间
	const auto merged_face_min_collision_patch = EvaluateMergedPatchToolCollision(
		slicer_load_merged_patch,
		merged_face_source_patch_id,
		cutting_tool);

	std::cout << "[Info] merged-face min-collision-patch matrix computed: faces="
		<< merged_face_min_collision_patch.size() << ", orientations="
		<< sampling_subtractive.sample_points.size() << std::endl;

	// ---------------- graph cut on merged patch ----------------
	{
		const int face_count = static_cast<int>(slicer_load_merged_patch.triangles.size());
		const int ori_count = static_cast<int>(sampling_subtractive.sample_points.size());
		const int block_count = graph_cut_block_count;
		const int num_labels = block_count * ori_count;
		const int INF_COST = 10000000;

		if (face_count == 0 || ori_count == 0 || block_count <= 0) {
			std::cout << "[Warn] skip graph cut: invalid sizes. face_count=" << face_count
				<< ", ori_count=" << ori_count << ", block_count=" << block_count << std::endl;
		}
		else {
			auto encode_label = [ori_count](int patch_id, int ori_id) -> int {
				// patch_id: 1..block_count, ori_id: 0..ori_count-1
				return (patch_id - 1) * ori_count + ori_id;
				};

			auto decode_label = [ori_count](int label_id, int& patch_id, int& ori_id) {
				patch_id = label_id / ori_count + 1;
				ori_id = label_id % ori_count;
				};

			// data cost: face x label
			std::vector<std::vector<int>> data_value(face_count, std::vector<int>(num_labels, INF_COST));

			int warn_a_lt_b = 0;
			int fallback_all_orientation_face_count = 0;

			for (int i = 0; i < face_count; ++i) {
				int b = 1;
				if (i < static_cast<int>(merged_face_source_patch_id.size())) {
					b = merged_face_source_patch_id[i];
				}
				else {
					std::cout << "[what?]i >= static_cast<int>(merged_face_source_patch_id.size())" << std::endl;
				}
				b = std::max(1, std::min(block_count, b));

				int feasible_cnt = 0;

				for (int ori = 0; ori < ori_count; ++ori) {
					int a = merged_face_min_collision_patch[i][ori];
					if (a == -1) {
						a = 1;
					}

					a = std::max(1, std::min(a, block_count + 1));

					// 按 [a..b] 添加标签，若 a>b 不添加并提示
					if (a > b) {
						++warn_a_lt_b;
						continue;
					}

					for (int p = a; p <= b; ++p) {
						const int lid = encode_label(p, ori);
						data_value[i][lid] = 0;
						++feasible_cnt;
					}
				}

				// 若碰撞约束没有产生可行 label，则固定为该面的来源 block，
				// 并开放其全部方向，让邻接和平滑项选择最终 orientation。
				if (feasible_cnt == 0) {
					++fallback_all_orientation_face_count;
					if (!EnableAllOrientationLabelsForSourceBlock(data_value[i], b, ori_count)) {
						std::cerr << "[Error] failed to enable source-block orientation labels: face="
							<< i << ", source_block=" << b << std::endl;
						return;
					}
				}
			}

			if (warn_a_lt_b > 0) {
				std::cout << "[Warn] graph-cut label range invalid (a>b) count = " << warn_a_lt_b << std::endl;
			}
			if (fallback_all_orientation_face_count > 0) {
				std::cout << "[Warn] faces with no collision-derived feasible label = "
					<< fallback_all_orientation_face_count
					<< " (all orientations of each face's source block enabled)." << std::endl;
			}

			MergedPatchGraphCutAdjacency graph_adjacency;
			if (!BuildMergedPatchGraphCutAdjacency(
				slicer_load_merged_patch,
				merged_face_source_patch_id,
				graph_cut_block_boundary_overlap_cost_per_unit,
				graph_adjacency)) {
				std::cerr << "[GraphCut] Skip global decomposition because merged-patch "
					<< "adjacency construction failed." << std::endl;
				return;
			}
			const auto& pixels_relations = graph_adjacency.face_relations;
			const auto& length_edges = graph_adjacency.edge_weights_times_ten;

			std::cout << "[Info] graph-cut nodes=" << face_count
				<< ", labels=" << num_labels
				<< ", edges=" << pixels_relations.size() << std::endl;

			std::vector<int> local_initial_labels;
			if (initialize_with_local_labels) {
				// Temporarily reduce each global feasible interval [a, b] to the
				// local greedy choice {a}. Restoring these entries to zero below
				// recovers the exact global data costs without another collision pass.
				for (int face_id = 0; face_id < face_count; ++face_id) {
					int source_block = 1;
					if (face_id < static_cast<int>(merged_face_source_patch_id.size())) {
						source_block = merged_face_source_patch_id[face_id];
					}
					source_block = std::max(1, std::min(block_count, source_block));

					for (int orientation_id = 0; orientation_id < ori_count; ++orientation_id) {
						int first_feasible_block =
							merged_face_min_collision_patch[face_id][orientation_id];
						if (first_feasible_block == -1) {
							first_feasible_block = 1;
						}
						first_feasible_block =
							std::max(1, std::min(first_feasible_block, block_count + 1));
						if (first_feasible_block > source_block) {
							continue;
						}

						for (int patch_id = first_feasible_block + 1;
							patch_id <= source_block;
							++patch_id) {
							data_value[face_id][encode_label(patch_id, orientation_id)] = INF_COST;
						}
					}
				}

				std::cout << "[GraphCut] Stage 1/2: solve local greedy labeling." << std::endl;
				local_initial_labels = GeneralGraph_DArraySArraySpatVarying(
					face_count,
					num_labels,
					data_value,
					pixels_relations,
					length_edges);
				auto local_label_types = local_initial_labels;
				std::sort(local_label_types.begin(), local_label_types.end());
				const int local_label_count = static_cast<int>(
					std::unique(local_label_types.begin(), local_label_types.end())
					- local_label_types.begin());
				std::cout << "[GraphCut] Local initialization uses "
					<< local_label_count << " labels." << std::endl;

				for (int face_id = 0; face_id < face_count; ++face_id) {
					int source_block = 1;
					if (face_id < static_cast<int>(merged_face_source_patch_id.size())) {
						source_block = merged_face_source_patch_id[face_id];
					}
					source_block = std::max(1, std::min(block_count, source_block));

					for (int orientation_id = 0; orientation_id < ori_count; ++orientation_id) {
						int first_feasible_block =
							merged_face_min_collision_patch[face_id][orientation_id];
						if (first_feasible_block == -1) {
							first_feasible_block = 1;
						}
						first_feasible_block =
							std::max(1, std::min(first_feasible_block, block_count + 1));
						if (first_feasible_block > source_block) {
							continue;
						}

						for (int patch_id = first_feasible_block + 1;
							patch_id <= source_block;
							++patch_id) {
							data_value[face_id][encode_label(patch_id, orientation_id)] = 0;
						}
					}
				}
			}

			std::vector<int> result_labels;
			if (initialize_with_local_labels) {
				std::cout << "[GraphCut] Stage 2/2: solve global labeling from the local result."
					<< std::endl;
				const auto initialized_result = RunGeneralGraphCutWithInitialLabels(
					face_count,
					num_labels,
					data_value,
					pixels_relations,
					length_edges,
					local_initial_labels,
					"local result");
				result_labels = initialized_result.labels;
			}
			else if (use_random_multistart) {
				const InitializedGraphCutResult multistart_result =
					RunRandomFeasibleLabelGraphCutRestarts(
						face_count,
						num_labels,
						data_value,
						pixels_relations,
						length_edges,
						random_restart_count,
						random_seed,
						random_thread_count,
						"global",
						PatchVisPath(
							"random_multistart_global_graphcut_runs.csv"));
				if (multistart_result.succeeded) {
					result_labels = multistart_result.labels;
				}
				else {
					std::cerr << "[GraphCut] All random restarts failed; use the default global initialization."
						<< std::endl;
					result_labels = GeneralGraph_DArraySArraySpatVarying(
						face_count,
						num_labels,
						data_value,
						pixels_relations,
						length_edges);
				}
			}
			else {
				result_labels = GeneralGraph_DArraySArraySpatVarying(
					face_count,
					num_labels,
					data_value,
					pixels_relations,
					length_edges);
			}
			auto final_label_types = result_labels;
			std::sort(final_label_types.begin(), final_label_types.end());
			const int final_label_count = static_cast<int>(
				std::unique(final_label_types.begin(), final_label_types.end())
				- final_label_types.begin());
			std::cout << "[GraphCut] Final "
				<< (initialize_with_local_labels
					? "local-initialized global"
					: (use_random_multistart ? "random-multistart global" : "global"))
				<< " result uses " << final_label_count << " labels." << std::endl;
			WriteMergedPatchGraphCutEdgeCostReport(
				PatchVisPath(
					output_prefix + "merged_patch_graphcut_edge_costs.csv"),
				PatchVisPath(
					output_prefix
					+ "merged_patch_graphcut_boundary_overlap_summary.txt"),
				graph_adjacency,
				result_labels);

			// 可视化：按最终label所属patch着色（离散固定色）
			const std::string gc_obj_file =
				PatchVisPath(output_prefix + "merged_patch_graphcut_label.obj");
			std::ofstream ofs(gc_obj_file);
			if (!ofs.is_open()) {
				std::cout << "[Warn] cannot open file for writing: " << gc_obj_file << std::endl;
			}
			else {
				auto color_from_patch = [](int ori_id) -> std::array<double, 3> {
					static const std::array<std::array<double, 3>, 32> palette = { {
						{0.894, 0.102, 0.110}, {0.216, 0.494, 0.722}, {0.302, 0.686, 0.290}, {0.596, 0.306, 0.639},
						{1.000, 0.498, 0.000}, {1.000, 1.000, 0.200}, {0.651, 0.337, 0.157}, {0.969, 0.506, 0.749},
						{0.600, 0.600, 0.600}, {0.121, 0.466, 0.705}, {0.682, 0.780, 0.909}, {1.000, 0.733, 0.470},
						{0.172, 0.627, 0.172}, {0.839, 0.153, 0.157}, {0.580, 0.404, 0.741}, {0.549, 0.337, 0.294},
						{0.890, 0.467, 0.761}, {0.498, 0.498, 0.498}, {0.737, 0.741, 0.133}, {0.090, 0.745, 0.811},
						{0.400, 0.760, 0.647}, {0.988, 0.553, 0.384}, {0.553, 0.627, 0.796}, {0.906, 0.541, 0.765},
						{0.651, 0.847, 0.329}, {1.000, 0.851, 0.184}, {0.898, 0.768, 0.580}, {0.702, 0.702, 0.702},
						{0.984, 0.603, 0.600}, {0.800, 0.922, 0.773}, {0.871, 0.796, 0.894}, {0.996, 0.851, 0.651}
					} };

					int idx = ori_id;
					if (idx < static_cast<int>(palette.size())) {
						return palette[static_cast<size_t>(idx)];
					}

					return palette[idx % palette.size()];

					};

				// 同时输出 face -> (patch,ori)
				const std::string gc_label_file =
					PatchVisPath(output_prefix + "merged_patch_graphcut_label.txt");
				std::ofstream lfs(gc_label_file);
				lfs << "face_id patch_id ori_id label_id\n";

				int v_count = 0;
				int skipped_faces = 0;
				int invalid_vertices = 0;

				// 按 patch_id 收集面，供 polyscope 分 mesh 可视化
				std::map<int, std::vector<vasco::core::Tri3>> patch_to_faces;
				std::map<int, std::vector<vasco::core::Tri3>> patch_id_to_faces;
				std::map<std::pair<int, int>, std::vector<vasco::core::Tri3>> patch_ori_to_faces;

				for (int i = 0; i < face_count; ++i) {
					int patch_id = 1, ori_id = 0;
					decode_label(result_labels[i], patch_id, ori_id);

					// 修正：按 patch_id 取颜色
					const auto color = color_from_patch(result_labels[i]);
					const auto& tri = slicer_load_merged_patch.triangles[i];

					// tri 索引与坐标合法性检查
					bool face_ok = true;
					for (int k = 0; k < 3; ++k) {
						const int vid = tri[k];
						if (vid < 0 || vid >= static_cast<int>(slicer_load_merged_patch.positions.size())) {
							face_ok = false;
							break;
						}
						const auto& p = slicer_load_merged_patch.positions[vid];
						if (!std::isfinite(p[0]) || !std::isfinite(p[1]) || !std::isfinite(p[2])) {
							face_ok = false;
							++invalid_vertices;
							break;
						}
					}

					if (!face_ok) {
						++skipped_faces;
						continue;
					}

					for (int k = 0; k < 3; ++k) {
						const auto& p = slicer_load_merged_patch.positions[tri[k]];
						ofs << "v " << p[0] << " " << p[1] << " " << p[2] << " "
							<< color[0] << " " << color[1] << " " << color[2] << "\n";
					}
					ofs << "f " << (v_count + 1) << " " << (v_count + 2) << " " << (v_count + 3) << "\n";
					v_count += 3;

					lfs << i << " " << patch_id << " " << ori_id << " " << result_labels[i] << "\n";

					// 保存到 patch 子网格
					patch_to_faces[result_labels[i]].push_back(tri);
					patch_id_to_faces[patch_id].push_back(tri);
					patch_ori_to_faces[{ patch_id, ori_id }].push_back(tri);
				}

				if (skipped_faces > 0) {
					std::cout << "[Warn] skipped faces while writing OBJ: " << skipped_faces
						<< ", invalid_vertices=" << invalid_vertices << std::endl;
				}

				for (const auto& kv : patch_id_to_faces) {
					const int patch_id = kv.first;
					const auto& patch_faces = kv.second;
					if (patch_faces.empty()) {
						continue;
					}

					const std::string patch_obj_file =
						PatchVisPath(
							output_prefix + "merged_patch_graphcut_patch_id_"
							+ std::to_string(patch_id)
							+ ".obj");
					std::ofstream pofs(patch_obj_file);
					if (!pofs.is_open()) {
						std::cout << "[Warn] cannot open file for writing: " << patch_obj_file << std::endl;
						continue;
					}

					const auto color = color_from_patch(patch_id);
					int patch_v_count = 0;
					for (const auto& tri : patch_faces) {
						for (int k = 0; k < 3; ++k) {
							const auto& p = slicer_load_merged_patch.positions[tri[k]];
							pofs << "v " << p[0] << " " << p[1] << " " << p[2] << " "
								<< color[0] << " " << color[1] << " " << color[2] << "\n";
						}
						pofs << "f " << (patch_v_count + 1) << " " << (patch_v_count + 2) << " " << (patch_v_count + 3) << "\n";
						patch_v_count += 3;
					}

					std::cout << "[Info] wrote graph-cut patch OBJ: " << patch_obj_file
						<< ", faces=" << patch_faces.size() << std::endl;
				}

				for (const auto& kv : patch_ori_to_faces) {
					const int patch_id = kv.first.first;
					const int ori_id = kv.first.second;
					const auto& patch_faces = kv.second;
					if (patch_faces.empty()) {
						continue;
					}

					const std::string patch_ori_obj_file =
						PatchVisPath(output_prefix + "merged_patch_graphcut_patch_id_" + std::to_string(patch_id)
							+ "_ori_id_" + std::to_string(ori_id) + ".obj");
					std::ofstream pofs(patch_ori_obj_file);
					if (!pofs.is_open()) {
						std::cout << "[Warn] cannot open file for writing: " << patch_ori_obj_file << std::endl;
						continue;
					}

					const auto color = color_from_patch(ori_id);
					int patch_v_count = 0;
					for (const auto& tri : patch_faces) {
						for (int k = 0; k < 3; ++k) {
							const auto& p = slicer_load_merged_patch.positions[tri[k]];
							pofs << "v " << p[0] << " " << p[1] << " " << p[2] << " "
								<< color[0] << " " << color[1] << " " << color[2] << "\n";
						}
						pofs << "f " << (patch_v_count + 1) << " " << (patch_v_count + 2) << " " << (patch_v_count + 3) << "\n";
						patch_v_count += 3;
					}

					std::cout << "[Info] wrote graph-cut patch/orientation OBJ: " << patch_ori_obj_file
						<< ", faces=" << patch_faces.size() << std::endl;

					Eigen::Vector3d tool_dir(0.0, 0.0, 1.0);
					if (ori_id >= 0 && ori_id < static_cast<int>(sampling_subtractive.sample_points.size())) {
						tool_dir = sampling_subtractive.sample_points[ori_id];
					}
					const double tool_dir_norm = tool_dir.norm();
					if (tool_dir_norm > 1e-12) {
						tool_dir /= tool_dir_norm;
					}

					const std::string tool_dir_file =
						PatchVisPath(output_prefix + "merged_patch_graphcut_patch_id_" + std::to_string(patch_id)
							+ "_ori_id_" + std::to_string(ori_id) + ".txt");
					std::ofstream tdfs(tool_dir_file);
					if (!tdfs.is_open()) {
						std::cout << "[Warn] cannot open file for writing: " << tool_dir_file << std::endl;
					}
					else {
						tdfs << tool_dir.x() << " " << tool_dir.y() << " " << tool_dir.z() << "\n";
						std::cout << "[Info] wrote graph-cut patch/orientation tool direction: "
							<< tool_dir_file << std::endl;
					}
				}

				// polyscope: 按 patch_id 拆分 mesh，并设置固定颜色
				std::vector<Eigen::Vector3d> patch_arrow_points;
				std::vector<Eigen::Vector3d> patch_arrow_dirs;
				std::vector<std::vector<double>> patch_arrow_colors;

				Eigen::Vector3d bb_min(MAX_D, MAX_D, MAX_D), bb_max(MIN_D, MIN_D, MIN_D);
				for (const auto& p : slicer_load_merged_patch.positions) {
					bb_min.x() = std::min(bb_min.x(), p[0]);
					bb_min.y() = std::min(bb_min.y(), p[1]);
					bb_min.z() = std::min(bb_min.z(), p[2]);
					bb_max.x() = std::max(bb_max.x(), p[0]);
					bb_max.y() = std::max(bb_max.y(), p[1]);
					bb_max.z() = std::max(bb_max.z(), p[2]);
				}
				const double model_diag = (bb_max - bb_min).norm();
				const double arrow_start_offset = std::max(1.0, 0.06 * model_diag);
				const double arrow_length = std::max(2.0, 0.12 * model_diag);

				for (const auto& kv : patch_to_faces) {
					const int patch_id = kv.first;
					const auto& patch_faces = kv.second;
					if (patch_faces.empty()) continue;
					int de_patch;
					int ori_id = 0;
					decode_label(patch_id, de_patch, ori_id);
					const std::string mesh_name =
						output_prefix + "merged_patch_graphcut_patch_" + std::to_string(height_of_beam_search) + "_" + std::to_string(de_patch) + "_" + std::to_string(ori_id);

					auto* ps_mesh = polyscope::registerSurfaceMesh(
						mesh_name,
						slicer_load_merged_patch.positions,
						patch_faces);

					const auto c = color_from_patch(ori_id);
					ps_mesh->setSurfaceColor({ (float)c[0], (float)c[1], (float)c[2] });

					std::cout << "[debug]mesh_name = " << mesh_name << ", patch_id=" << patch_id << ", de_patch=" << de_patch << ", ori_id=" << ori_id
						<< ", c= " << c[0] << "," << c[1] << "," << c[2] << std::endl;

					Eigen::Vector3d patch_center(0.0, 0.0, 0.0);
					int center_count = 0;
					for (const auto& tri : patch_faces) {
						for (int k = 0; k < 3; ++k) {
							const int vid = tri[k];
							if (vid >= 0 && vid < static_cast<int>(slicer_load_merged_patch.positions.size())) {
								const auto& p = slicer_load_merged_patch.positions[vid];
								patch_center.x() += p[0];
								patch_center.y() += p[1];
								patch_center.z() += p[2];
								++center_count;
							}
						}
					}
					if (center_count > 0) {
						patch_center /= static_cast<double>(center_count);

						Eigen::Vector3d ori_dir(0.0, 0.0, 1.0);
						if (ori_id >= 0 && ori_id < static_cast<int>(sampling_subtractive.sample_points.size())) {
							ori_dir = sampling_subtractive.sample_points[ori_id];
						}
						Eigen::Vector3d ori_vec(ori_dir);
						double n = ori_vec.norm();
						if (n < 1e-12) {
							ori_vec = Eigen::Vector3d(0.0, 0.0, 1.0);
						}
						else {
							ori_vec /= n;
						}

						const Eigen::Vector3d arrow_start = patch_center + ori_vec * arrow_start_offset;
						const Eigen::Vector3d arrow_end = arrow_start + ori_vec * arrow_length;

						patch_arrow_points.push_back(arrow_start);
						patch_arrow_dirs.push_back(ori_vec);

						patch_arrow_colors.push_back({ c[0], c[1], c[2] });

						const std::string arrow_name =
							output_prefix + "merged_patch_graphcut_arrow_" + std::to_string(height_of_beam_search) + "_" + std::to_string(de_patch) + "_" + std::to_string(ori_id);

						Eigen::Vector3d ref_axis = (std::abs(ori_vec.z()) < 0.9)
							? Eigen::Vector3d(0.0, 0.0, 1.0)
							: Eigen::Vector3d(0.0, 1.0, 0.0);
						Eigen::Vector3d perp1 = ori_vec.cross(ref_axis);
						double perp1_norm = perp1.norm();
						if (perp1_norm < 1e-12) {
							perp1 = Eigen::Vector3d(1.0, 0.0, 0.0);
						}
						else {
							perp1 /= perp1_norm;
						}

						const double shaft_w = std::max(0.2, arrow_length * 0.08);
						const double head_len = arrow_length * 0.35;
						const double head_w = shaft_w * 2.0;

						Eigen::Vector3d m = arrow_end - ori_vec * head_len;
						Eigen::Vector3d s0 = arrow_start + perp1 * shaft_w;
						Eigen::Vector3d s1 = arrow_start - perp1 * shaft_w;
						Eigen::Vector3d m0 = m + perp1 * shaft_w;
						Eigen::Vector3d m1 = m - perp1 * shaft_w;
						Eigen::Vector3d h0 = m + perp1 * head_w;
						Eigen::Vector3d h1 = m - perp1 * head_w;

						std::vector<vasco::core::Vec3> arrow_pos = {
							{ s0.x(), s0.y(), s0.z() },
							{ s1.x(), s1.y(), s1.z() },
							{ m1.x(), m1.y(), m1.z() },
							{ m0.x(), m0.y(), m0.z() },
							{ h0.x(), h0.y(), h0.z() },
							{ h1.x(), h1.y(), h1.z() },
							{ arrow_end.x(), arrow_end.y(), arrow_end.z() }
						};
						std::vector<vasco::core::Tri3> arrow_tri = {
							{ 0, 1, 2 }, { 0, 2, 3 },
							{ 4, 5, 6 }, { 3, 4, 6 }, { 2, 6, 5 }
						};

						auto* ps_arrow = polyscope::registerSurfaceMesh(arrow_name, arrow_pos, arrow_tri);
						ps_arrow->setSurfaceColor({ (float)c[0], (float)c[1], (float)c[2] });
					}
				}

				if (!patch_arrow_points.empty()) {
					Visual vis_patch_ori;
					vis_patch_ori.generateModelForRendering_11(
						patch_arrow_points,
						patch_arrow_dirs,
						patch_arrow_colors,
						PatchVisPath(output_prefix + "merged_patch_graphcut_patch_ori_arrows.obj"));
				}
				polyscope::show();

				std::cout << "[Info] wrote graph-cut colored OBJ: " << gc_obj_file << std::endl;
				std::cout << "[Info] wrote graph-cut labels txt: " << gc_label_file << std::endl;
			}
		}
	}

}

int HybridManufacturing::subtractive_accessibility_decomposition_local(int height_of_beam_search, cutter cutting_tool)
{
	return RunSubtractiveAccessibilityDecompositionLocal(
		height_of_beam_search,
		cutting_tool,
		0,
		0,
		1);
}

int HybridManufacturing::subtractive_accessibility_decomposition_random_multistart_local(
	int height_of_beam_search,
	cutter cutting_tool,
	int restart_count,
	unsigned int random_seed,
	int thread_count)
{
	return RunSubtractiveAccessibilityDecompositionLocal(
		height_of_beam_search,
		cutting_tool,
		std::max(1, restart_count),
		random_seed,
		thread_count);
}

int HybridManufacturing::RunSubtractiveAccessibilityDecompositionLocal(
	int height_of_beam_search,
	cutter cutting_tool,
	int random_restart_count,
	unsigned int random_seed,
	int random_thread_count)
{
	const bool use_random_multistart = random_restart_count > 0;
	const std::string output_prefix =
		use_random_multistart ? "random_multistart_local_" : "";

	sampling_subtractive.OrientationSamplePoints();	//sampling_subtractive生成球面采样点
	cutting_tool.cylinder_r = 1.5;
	cutting_tool.cylinder_height = 28;
	cutting_tool.ball_r = 1.5;
	cutting_tool.carriage_r = 23;
	cutting_tool.carriage_height = 33;

	/////show all accessible points in every orientation/////
	Eigen::MatrixXd V_B2;
	Eigen::MatrixXi F_B2;
	igl::readOBJ("ball.obj", V_B2, F_B2);

	Slicer_2 slicer_load_current_patch;

	Slicer_2 slicer_load_next_patch;

	const std::string ancestor_source_report_file =
		AncestorSourceVisPath(
			"beam_search_last_layer_ancestor_sources_"
			+ MakeAccessibilityDebugFileToken(file_name)
			+ ".txt");
	const std::vector<std::string> ancestor_patch_files =
		LoadFirstFinalNodeAncestorPatchFiles(
			ancestor_source_report_file,
			file_name,
			PatchVisDir());
	const int graph_cut_block_count = ancestor_patch_files.empty()
		? height_of_beam_search
		: static_cast<int>(ancestor_patch_files.size());

	if (!ancestor_patch_files.empty()) {
		std::cout << "[Info] local decomposition uses ancestor patches from first final node in "
			<< ancestor_source_report_file << std::endl;
		for (int i = 0; i < static_cast<int>(ancestor_patch_files.size()); ++i) {
			std::cout << "  ancestor_patch[" << (i + 1) << "]: "
				<< ancestor_patch_files[i] << std::endl;
		}
	}

	// 与 global 使用同一条祖先 patch 链；报告不存在时使用相同的旧格式回退。
	std::vector<int> merged_vertex_source_patch_id;
	std::vector<int> merged_face_source_patch_id;

	Slicer_2 slicer_load_merged_patch = ancestor_patch_files.empty()
		? MergeBlockPatchesWithDedup(
			height_of_beam_search,
			merged_vertex_source_patch_id,
			merged_face_source_patch_id,
			1e-3)
		: MergeBlockPatchesWithDedup(
			ancestor_patch_files,
			merged_vertex_source_patch_id,
			merged_face_source_patch_id,
			1e-3);
	if (!StitchMergedPatchBoundariesForGraphCut(
		slicer_load_merged_patch,
		merged_face_source_patch_id,
		use_random_multistart ? "random-multistart local" : "local")) {
		return 999;
	}
	slicer_load_merged_patch.save(
		PatchVisPath(
			use_random_multistart
			? "random_multistart_local_block_patch_merged_removedup.obj"
			: "block_patch_merged_removedup_local.obj"));
	std::cout << "[Info] merged patch mesh: V=" << slicer_load_merged_patch.positions.size()
		<< ", F=" << slicer_load_merged_patch.triangles.size() << std::endl;

	// 计算 merged patch 上每个三角面在每个采样方向下的最大碰撞 patch_index
	const auto merged_face_min_collision_patch = EvaluateMergedPatchToolCollision(
		slicer_load_merged_patch,
		merged_face_source_patch_id,
		cutting_tool);

	std::cout << "[Info] merged-face min-collision-patch matrix computed: faces="
		<< merged_face_min_collision_patch.size() << ", orientations="
		<< sampling_subtractive.sample_points.size() << std::endl;

	int ret = 999;
	// ---------------- graph cut on merged patch ----------------
	{
		const int face_count = static_cast<int>(slicer_load_merged_patch.triangles.size());
		const int ori_count = static_cast<int>(sampling_subtractive.sample_points.size());
		const int block_count = graph_cut_block_count;
		const int num_labels = block_count * ori_count;
		const int INF_COST = 10000000;

		if (face_count == 0 || ori_count == 0 || block_count <= 0) {
			std::cout << "[Warn] skip graph cut: invalid sizes. face_count=" << face_count
				<< ", ori_count=" << ori_count << ", block_count=" << block_count << std::endl;
		}
		else {
			auto encode_label = [ori_count](int patch_id, int ori_id) -> int {
				// patch_id: 1..block_count, ori_id: 0..ori_count-1
				return (patch_id - 1) * ori_count + ori_id;
				};

			// data cost: face x label
			std::vector<std::vector<int>> data_value(face_count, std::vector<int>(num_labels, INF_COST));

			int warn_a_lt_b = 0;
			int fallback_all_orientation_face_count = 0;

			for (int i = 0; i < face_count; ++i) {
				int b = 1;
				if (i < static_cast<int>(merged_face_source_patch_id.size())) {
					b = merged_face_source_patch_id[i];
				}
				else {
					std::cout << "[what?]i >= static_cast<int>(merged_face_source_patch_id.size())" << std::endl;
				}
				b = std::max(1, std::min(block_count, b));

				int feasible_cnt = 0;

				for (int ori = 0; ori < ori_count; ++ori) {
					int a = merged_face_min_collision_patch[i][ori];
					if (a == -1) {
						a = 1;
					}

					a = std::max(1, std::min(a, block_count + 1));

					// Patch 编号按 beam-search 切除顺序递增，而实际制造时按相反顺序恢复 block。
					// 因此在可行区间 [a, b] 中只开放 a，表示尽量推迟该面在实际制造流程中的减材时刻。
					if (a > b) {
						++warn_a_lt_b;
						continue;
					}

					for (int p = a; p <= a; ++p) {
						const int lid = encode_label(p, ori);
						data_value[i][lid] = 0;
						++feasible_cnt;
					}
				}

				// 若碰撞约束没有产生可行 label，则固定为该面的来源 block，
				// 并开放其全部方向，让邻接和平滑项选择最终 orientation。
				if (feasible_cnt == 0) {
					++fallback_all_orientation_face_count;
					if (!EnableAllOrientationLabelsForSourceBlock(data_value[i], b, ori_count)) {
						std::cerr << "[Error] failed to enable source-block orientation labels: face="
							<< i << ", source_block=" << b << std::endl;
						return 999;
					}
				}
			}

			if (warn_a_lt_b > 0) {
				std::cout << "[Warn] graph-cut label range invalid (a>b) count = " << warn_a_lt_b << std::endl;
			}
			if (fallback_all_orientation_face_count > 0) {
				std::cout << "[Warn] faces with no collision-derived feasible label = "
					<< fallback_all_orientation_face_count
					<< " (all orientations of each face's source block enabled)." << std::endl;
			}

			MergedPatchGraphCutAdjacency graph_adjacency;
			if (!BuildMergedPatchGraphCutAdjacency(
				slicer_load_merged_patch,
				merged_face_source_patch_id,
				graph_cut_block_boundary_overlap_cost_per_unit,
				graph_adjacency)) {
				std::cerr << "[GraphCut] Skip local decomposition because merged-patch "
					<< "adjacency construction failed." << std::endl;
				return 999;
			}
			const auto& pixels_relations = graph_adjacency.face_relations;
			const auto& length_edges = graph_adjacency.edge_weights_times_ten;

			std::cout << "[Info] graph-cut nodes=" << face_count
				<< ", labels=" << num_labels
				<< ", edges=" << pixels_relations.size() << std::endl;

			// graph cut
			std::vector<int> result_labels;
			if (use_random_multistart) {
				const InitializedGraphCutResult multistart_result =
					RunRandomFeasibleLabelGraphCutRestarts(
						face_count,
						num_labels,
						data_value,
						pixels_relations,
						length_edges,
						random_restart_count,
						random_seed,
						random_thread_count,
						"local",
						PatchVisPath(
							output_prefix + "graphcut_runs.csv"));
				if (multistart_result.succeeded) {
					result_labels = multistart_result.labels;
				}
				else {
					std::cerr << "[GraphCut] All local random restarts failed; "
						<< "use the default local initialization." << std::endl;
					result_labels = GeneralGraph_DArraySArraySpatVarying(
						face_count,
						num_labels,
						data_value,
						pixels_relations,
						length_edges);
				}
			}
			else {
				result_labels = GeneralGraph_DArraySArraySpatVarying(
					face_count,
					num_labels,
					data_value,
					pixels_relations,
					length_edges);
			}

			auto label_type_count = result_labels;
			std::sort(label_type_count.begin(), label_type_count.end());
			ret = std::unique(label_type_count.begin(), label_type_count.end()) - label_type_count.begin();

			const std::string visualization_prefix =
				use_random_multistart ? output_prefix : "local_";
			WriteMergedPatchGraphCutEdgeCostReport(
				PatchVisPath(
					visualization_prefix + "merged_patch_graphcut_edge_costs.csv"),
				PatchVisPath(
					visualization_prefix
					+ "merged_patch_graphcut_boundary_overlap_summary.txt"),
				graph_adjacency,
				result_labels);
			accessibility_visualization::ShowMergedPatchGraphCutLabels(
				PatchVisDir(),
				visualization_prefix,
				height_of_beam_search,
				ori_count,
				slicer_load_merged_patch,
				result_labels,
				sampling_subtractive.sample_points);
		}
	}
	return ret;
}

vector<vector<int>> HybridManufacturing::getAccessOri(const Slicer_2& slicer, Slicer_2& slicer_load_patch, vector<vasco::core::Vec3>& all_sample_points_in_triangles, cutter cutting_tool)
{
	auto& need_detect_triangle = slicer_load_patch.triangles;
	vector<vector<int>> accessible_ori_of_need_detect_V(need_detect_triangle.size(), vector<int>(sampling_subtractive.sample_points.size(), 10000000));
	cout << "%sampling_subtractive.sample_points.size() " << sampling_subtractive.sample_points.size() << endl;
	for (int ori = 0; ori < sampling_subtractive.sample_points.size(); ori++) {
		vector<std::vector<Eigen::Vector3d>> temp_new_V_remain(slicer.triangles.size());
		for (int i = 0; i < slicer.triangles.size(); i++) {
			temp_new_V_remain[i].resize(3);
			for (int k = 0; k < 3; k++) {
				temp_new_V_remain[i][k](0, 0) = slicer.positions[slicer.triangles[i][k]][0];
				temp_new_V_remain[i][k](1, 0) = slicer.positions[slicer.triangles[i][k]][1];
				temp_new_V_remain[i][k](2, 0) = slicer.positions[slicer.triangles[i][k]][2];
			}
		}

		std::vector<Eigen::Vector3d> temp_V_need_detect(need_detect_triangle.size());
		for (int i = 0; i < need_detect_triangle.size(); i++) {
			temp_V_need_detect[i](0, 0) = all_sample_points_in_triangles[i][0];
			temp_V_need_detect[i](1, 0) = all_sample_points_in_triangles[i][1];
			temp_V_need_detect[i](2, 0) = all_sample_points_in_triangles[i][2];
		}

		vector<std::vector<Eigen::Vector3d>> temp_V_need_detect_triangle(need_detect_triangle.size());
		for (int i = 0; i < need_detect_triangle.size(); i++) {
			temp_V_need_detect_triangle[i].resize(3);
			for (int k = 0; k < 3; k++) {
				temp_V_need_detect_triangle[i][k](0, 0) = slicer_load_patch.positions[need_detect_triangle[i][k]][0];
				temp_V_need_detect_triangle[i][k](1, 0) = slicer_load_patch.positions[need_detect_triangle[i][k]][1];
				temp_V_need_detect_triangle[i][k](2, 0) = slicer_load_patch.positions[need_detect_triangle[i][k]][2];
			}
		}

		Eigen::Matrix3d rotMatrix;
		Eigen::Vector3d vectorBefore(0, 0, 1);
		Eigen::Vector3d vectorAfter(sampling_subtractive.sample_points[ori]);
		vectorAfter.normalize();
		rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vectorAfter).toRotationMatrix();
		for (int i = 0; i < temp_new_V_remain.size(); i++)
			for (int j = 0; j < temp_new_V_remain[i].size(); j++)
				temp_new_V_remain[i][j] = rotMatrix.inverse() * temp_new_V_remain[i][j];
		for (int i = 0; i < temp_V_need_detect.size(); i++)
			temp_V_need_detect[i] = rotMatrix.inverse() * temp_V_need_detect[i];
		for (int i = 0; i < temp_V_need_detect_triangle.size(); i++)
			for (int j = 0; j < temp_V_need_detect_triangle[i].size(); j++)
				temp_V_need_detect_triangle[i][j] = rotMatrix.inverse() * temp_V_need_detect_triangle[i][j];

		vector<Eigen::Vector3d> all_normals_of_need_detect_triangle;
		all_normals_of_need_detect_triangle.clear();
		for (int i = 0; i < need_detect_triangle.size(); i++) {
			cv::Point3d Point1(temp_V_need_detect_triangle[i][0](0, 0), temp_V_need_detect_triangle[i][0](1, 0), temp_V_need_detect_triangle[i][0](2, 0));
			cv::Point3d Point2(temp_V_need_detect_triangle[i][1](0, 0), temp_V_need_detect_triangle[i][1](1, 0), temp_V_need_detect_triangle[i][1](2, 0));
			cv::Point3d Point3(temp_V_need_detect_triangle[i][2](0, 0), temp_V_need_detect_triangle[i][2](1, 0), temp_V_need_detect_triangle[i][2](2, 0));
			double na = (Point2.y - Point1.y) * (Point3.z - Point1.z) - (Point2.z - Point1.z) * (Point3.y - Point1.y);
			double nb = (Point2.z - Point1.z) * (Point3.x - Point1.x) - (Point2.x - Point1.x) * (Point3.z - Point1.z);
			double nc = (Point2.x - Point1.x) * (Point3.y - Point1.y) - (Point2.y - Point1.y) * (Point3.x - Point1.x);
			Eigen::Vector3d vn(na, nb, nc);
			vn.normalize();
			all_normals_of_need_detect_triangle.push_back(vn);
			//cout << vn.x() << " "<<vn.y()<<" "<<vn.z() << endl;
		}

		vector<double> max_z_of_triangles(temp_new_V_remain.size());
		for (int i = 0; i < temp_new_V_remain.size(); i++) {
			max_z_of_triangles[i] = MIN_D;
			for (int j = 0; j < 3; j++)
				max_z_of_triangles[i] = max(max_z_of_triangles[i], temp_new_V_remain[i][j](2, 0));
		}

		///////////////////collision detection////////////////////////
		PrepareToolForCollision(cutting_tool);

		for (int i = 0; i < temp_V_need_detect.size(); i++) {
			Eigen::Vector3d center_point;

			center_point.x() = temp_V_need_detect[i](0, 0) + (cutting_tool.cylinder_r) * all_normals_of_need_detect_triangle[i].x();
			center_point.y() = temp_V_need_detect[i](1, 0) + (cutting_tool.cylinder_r) * all_normals_of_need_detect_triangle[i].y();
			center_point.z() = temp_V_need_detect[i](2, 0) + (cutting_tool.cylinder_r) * all_normals_of_need_detect_triangle[i].z();

			bool jud_collision = false;

			for (int ii = 0; ii < temp_new_V_remain.size(); ii++) {
				if (CheckToolCollisionWithCell(center_point, temp_new_V_remain[ii], max_z_of_triangles[ii], cutting_tool, 30.0, 3.0)) {
					jud_collision = true;
					break;
				}
			}

			if (jud_collision == false) {
				//vis_points.push_back(temp_V_need_detect[i]);
				accessible_ori_of_need_detect_V[i][ori] = 0;
			}
			//else
				//cout << "% " << endl;
		}
	}
	return accessible_ori_of_need_detect_V;
}

void HybridManufacturing::EvaluateCandidateNode(
	int& i,
	vector<int>& candidate_nodes,
	OrientationScores& all_calculated_value,
	vector<BeamTreeEntry>& tree_entries,
	const vector<vector<Eigen::Vector3d>>& save_ori,
	int height_of_beam_search,
	int cont_number_of_queue,
	const string& file_name,
	int now_last_node,
	nozzle the_nozzle)
{
	printf("\r[%d%%]>", i * 100 / std::max(1, static_cast<int>(candidate_nodes.size()) - 1));
	for (int j = 1; j <= i * 20 / candidate_nodes.size(); j++)
		cout << "▉";

	int index_of_pre_node = tree_entries[candidate_nodes[i]].pre_queue_index;
	vector<vector<cv::Point3d>> all_cut_layers;
	vector<int> all_cut_layers_dependency_layer;
	all_cut_layers.clear();
	all_cut_layers_dependency_layer.clear();
	for (int j = 0; j < tree_entries[candidate_nodes[i]].cut_layers.size(); j++) {
		int index_of_layers = tree_entries[candidate_nodes[i]].cut_layers[j];
		all_cut_layers.push_back(tree_entries[candidate_nodes[i]].layers[index_of_layers]);
		all_cut_layers_dependency_layer.push_back(tree_entries[candidate_nodes[i]].cut_layers_dependency_layer[j]);
	}

	Slicer_2 slicer_G;
	if (tree_entries[now_last_node].judge_continue == false)
		slicer_G.load(file_name + "-" + to_string(height_of_beam_search - 1) + "_" + to_string(index_of_pre_node) + ".obj");
	else
		slicer_G.load(file_name + "-" + to_string(height_of_beam_search - 1) + "_" + to_string(index_of_pre_node) + "_" + to_string(tree_entries[candidate_nodes[i]].continue_id - 1) + "_subblock.obj");

	all_calculated_value[i] = GainMesh(
		slicer_G,
		all_cut_layers,
		save_ori[candidate_nodes[i]][0],
		height_of_beam_search,
		cont_number_of_queue,
		index_of_pre_node,
		all_cut_layers_dependency_layer,
		tree_entries[now_last_node].judge_continue,
		tree_entries[candidate_nodes[i]].continue_id);

	calculate_fragile_value(all_calculated_value[i], all_cut_layers, tree_entries[candidate_nodes[i]].fragile_v);

	tree_entries[candidate_nodes[i]].larger_base = all_calculated_value[i].large_base; //有点疑问，看看candidate_nodes[i]是否一定等于i
	if (all_calculated_value[i].value_of_self_support == 0) {
		candidate_nodes.erase(candidate_nodes.begin() + i);
		all_calculated_value.erase(all_calculated_value.begin() + i);
		i--;
		return;
	}

	detect_collision_with_printing_platform(
		i,
		candidate_nodes,
		all_calculated_value,
		all_cut_layers,
		save_ori[candidate_nodes[i]][0],
		the_nozzle);
}

void HybridManufacturing::LogSelectedCandidateMetrics(
	int cont_w,
	const vector<int>& candidate_nodes,
	int now_last_node,
	const vector<vector<Eigen::Vector3d>>& save_ori,
	const vector<BeamTreeEntry>& tree_entries,
	queue<int>& last_step_nodes,
	const OrientationScores& all_calculated_value,
	const OrientationScores& pure_value,
	int& sum_candidate_blocks,
	int& sum_connected_components,
	int& cont_extra_additive_orientation,
	vector<double>& evaluation_value)
{
	cout << endl << "selected orientation: " << save_ori[candidate_nodes[cont_w]][0].x() << " " << save_ori[candidate_nodes[cont_w]][0].y() << " " << save_ori[candidate_nodes[cont_w]][0].z();
	cout << endl << "number of candidate_nodes: " << candidate_nodes.size() << endl;
	sum_candidate_blocks += candidate_nodes.size();
	sum_connected_components += tree_entries[candidate_nodes[cont_w]].cut_layers.size();
	if (tree_entries[now_last_node].judge_continue == false)
		last_step_nodes.push(candidate_nodes[cont_w]);
	cout << "self support value of selected node: " << all_calculated_value[cont_w].value_of_self_support << endl << "■■■■■■■■■■■■■■■■■■■" << endl;
	cout << "value_of_more_slice_layers: " << all_calculated_value[cont_w].value_of_more_slice_layers << endl;
	cout << "value_of_area_S: " << all_calculated_value[cont_w].value_of_area_S << endl;
	cout << "value_of_covering_points: " << all_calculated_value[cont_w].value_of_covering_points << endl;
	cout << "value_of_less_clipping_plane: " << all_calculated_value[cont_w].value_of_less_clipping_plane << endl;
	cout << "value_of_fragile: " << all_calculated_value[cont_w].value_of_fragile << endl;
	cout << "value_of_orientation: " << all_calculated_value[cont_w].value_of_orientation << endl;
	cout << "value_of_projected: " << all_calculated_value[cont_w].value_of_projected << endl;
	cout << "■■■■■■■■■■■■■■■■■■■" << endl;

	cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
	cout << "Pure value orien:" << pure_value[cont_w].value_of_orientation << endl;
	cout << "Pure value fragile:" << pure_value[cont_w].value_of_fragile << endl;
	cout << "Pure value projection:" << pure_value[cont_w].value_of_projected << endl;
	cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
	if (tree_entries[candidate_nodes[cont_w]].judge_continue == true) {
		cout << "*****NEED CHANGE ORIENTATION*****" << endl;
		cont_extra_additive_orientation++;
	}

	evaluation_value[0] += all_calculated_value[cont_w].value_of_more_slice_layers;
	evaluation_value[1] += all_calculated_value[cont_w].value_of_covering_points;
	evaluation_value[2] += all_calculated_value[cont_w].value_of_less_clipping_plane;
	evaluation_value[3] += all_calculated_value[cont_w].value_of_orientation;
	evaluation_value[4] += all_calculated_value[cont_w].value_of_fragile;
	evaluation_value[5] += all_calculated_value[cont_w].value_of_projected;
}

void HybridManufacturing::BuildCutLayersForCandidate(
	int cont_w,
	const vector<int>& candidate_nodes,
	const vector<BeamTreeEntry>& tree_entries,
	vector<vector<cv::Point3d>>& all_cut_layers,
	vector<int>& all_cut_layers_dependency_layer,
	vector<int>& flag_cut_layers_is_hole) const
{
	all_cut_layers.clear();
	all_cut_layers_dependency_layer.clear();
	flag_cut_layers_is_hole.clear();
	const BeamTreeEntry& entry = tree_entries[candidate_nodes[cont_w]];
	for (int j = 0; j < static_cast<int>(entry.cut_layers.size()); ++j) {
		const int index_of_layers = entry.cut_layers[j];
		if (index_of_layers < 0
			|| index_of_layers >= static_cast<int>(entry.layers.size())
			|| entry.layers[index_of_layers].size() < 3) {
			continue;
		}
		const int dependency =
			j < static_cast<int>(entry.cut_layers_dependency_layer.size())
			? entry.cut_layers_dependency_layer[j]
			: 0;
		const int outer_cut_layer_id = static_cast<int>(all_cut_layers.size());
		all_cut_layers.push_back(entry.layers[index_of_layers]);
		all_cut_layers_dependency_layer.push_back(dependency);
		flag_cut_layers_is_hole.push_back(-1);

		if (index_of_layers < static_cast<int>(entry.holes.size())
			&& !entry.holes[index_of_layers].empty()) {
			for (const CutLayer& hole : entry.holes[index_of_layers]) {
				if (hole.size() < 3) {
					continue;
				}
				all_cut_layers.push_back(hole);
				all_cut_layers_dependency_layer.push_back(dependency);
				flag_cut_layers_is_hole.push_back(outer_cut_layer_id);
			}
		}
		else if (index_of_layers < static_cast<int>(entry.contain.size())
			&& entry.contain[index_of_layers].size() >= 3) {
			all_cut_layers.push_back(entry.contain[index_of_layers]);
			all_cut_layers_dependency_layer.push_back(dependency);
			flag_cut_layers_is_hole.push_back(outer_cut_layer_id);
		}
	}
}

bool HybridManufacturing::PrepareOuterBeamSearchNode(
	queue<int>& last_step_nodes,
	vector<vector<bool>>& is_fragile_V_2,
	int& now_last_node,
	int height_of_beam_search,
	int cont_number_of_queue,
	const vector<BeamTreeEntry>& tree_entries,
	double& sum_time_5)
{
	clock_t start_time_2 = clock();
	is_fragile_V_2 = is_fragile_V;	//fragile信息的快照
	now_last_node = last_step_nodes.front();
	/*if (Tree_nodes_error[now_last_node] == true) {
		last_step_nodes.pop();
		continue;
	}*/

	cout << endl << "Slicing and inner-beam-search......" << endl;
	//load blocks//
	const char* config_path = "config.ini";
	Katana::Instance().config.loadConfig(config_path);
	current_node_mesh.clear();
	if (tree_entries[now_last_node].judge_continue == false) {
		std::string current_file_name = file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + ".stl";
		Katana::Instance().stl.loadStl(current_file_name.c_str());	//加载当前节点对应的模型
		if (!vasco::mesh_validation::LoadAndCheckCurrentNodeMesh(current_file_name, current_node_mesh)) {	//加载当前节点对应的模型
			std::cout << "[BeamSearch] Skip node because LoadAndCheckCurrentNodeMesh failed: "
				<< "node=" << now_last_node
				<< ", source=" << current_file_name
				<< std::endl;
			clock_t end_time_2 = clock();
			sum_time_5 += double(end_time_2 - start_time_2) / CLOCKS_PER_SEC;
			return false;
		}

		cout << "T" << endl;
	}
	else {
		std::string current_file_name = file_name + "-" + to_string(height_of_beam_search - 1) + "_" + to_string(tree_entries[now_last_node].pre_queue_index) + ".stl";
		Katana::Instance().stl.loadStl(current_file_name.c_str());
		if (!vasco::mesh_validation::LoadAndCheckCurrentNodeMesh(current_file_name, current_node_mesh)) {
			std::cout << "[BeamSearch] Skip continue node because LoadAndCheckCurrentNodeMesh failed: "
				<< "node=" << now_last_node
				<< ", source=" << current_file_name
				<< std::endl;
			clock_t end_time_2 = clock();
			sum_time_5 += double(end_time_2 - start_time_2) / CLOCKS_PER_SEC;
			return false;
		}
		cout << "U" << now_last_node << endl;
	}
	Katana::Instance().temp_vertices.resize(Katana::Instance().vertices.size());	//更新当前模型的顶点信息
	Katana::Instance().temp_vertices = Katana::Instance().vertices;
	clock_t end_time_2 = clock();
	sum_time_5 += double(end_time_2 - start_time_2) / CLOCKS_PER_SEC;
	return true;
}

void HybridManufacturing::outer_beam_search(nozzle the_nozzle, cutter cutting_tool)
{
	int W1 = 4;  // 4
	clock_t start_time_total, end_time_total;
	clock_t start_time, end_time;
	clock_t start_time_6, end_time_6;
	clock_t start_time_7, end_time_7;
	clock_t start_time_8, end_time_8;
	float total_time_1 = 0, total_time_2 = 0, total_time_3 = 0, total_time_4 = 0;
	double sum_time = 0;
	double sum_time_2 = 0;
	double sum_time_3 = 0;
	double sum_time_4 = 0;
	double sum_time_5 = 0;
	double sum_time_6 = 0;
	int Sum_candidate_blocks = 0;
	int Sum_connected_components = 0;
	vector<double> evaluation_value(6, 0);
	cont_extra_additive_orientation = 0;

	start_time_total = clock();
	bool jud_outer_beam_search_terminate = false;
	vector<BeamTreeEntry> tree_entries;
	vector<int> candidate_nodes;	//当前层生成的候选节点
	queue<int> last_step_nodes;
	vector<vector<Eigen::Vector3d>> save_ori;
	const vector<int> root_cut_layers;
	const vector<Eigen::Vector3d> root_orientations{ Eigen::Vector3d::Zero() };
	BeamTreeEntry root_entry{};
	root_entry.parent_id = -1;
	root_entry.pre_queue_index = -1;
	root_entry.source_input_file = file_name + "-0_0.obj";
	tree_entries.push_back(root_entry);
	last_step_nodes.push(0);	//推入根节点，作为第一层搜索的初始节点
	tree_entries[0].judge_continue = false;
	tree_entries[0].continue_id = -1;
	tree_entries[0].error = false;
	save_ori.push_back(root_orientations);

	vector<vector<area_S>> ori_all_the_area_S = all_the_area_S;	//不可达点的all_the_area_S的拷贝
	ori_all_the_covering_points = all_the_covering_points;	//不可达点的all_the_covering_points的拷贝
	SAMPLE_ON_BALL sampling;
	int index_node = 0;
	sampling.OrientationSamplePoints_2();	//生成球面采样方向，作为增材打印方向候选
	//Visual Vis_ori;
	//Vis_ori.generateModelForRendering_6(sampling.sample_points, file_name);


	GetALLFragileVertex(sampling);

	Eigen::Matrix3d rotMatrix;
	vector<vector<int>> final_pathes_include_S;
	vector<vector<int>> final_pathes_include_sample_points;
	vector<vector<int>> final_pathes_include_covering_points;
	final_pathes_include_S.push_back(root_cut_layers);
	final_pathes_include_sample_points.push_back(root_cut_layers);
	final_pathes_include_covering_points.push_back(root_cut_layers);
	all_saved_mesh.resize(W1);

	int cont_number_of_queue = 0;	//记录当前处理到第几个节点，用于设定读取的文件名等
	int height_of_beam_search = 0;
	vector<int> last_selected_nodes;
	auto make_cut_source_file = [this](bool is_continue_block, int height, int pre_queue_index, int continue_id) {
		if (!is_continue_block) {
			return file_name + "-" + to_string(height - 1) + "_" + to_string(pre_queue_index) + ".obj";
		}
		return file_name + "-" + to_string(height - 1) + "_" + to_string(pre_queue_index)
			+ "_" + to_string(continue_id - 1) + "_subblock.obj";
		};

	ofstream ofile_cont("E:\\Hybrid manufacturing\\HybridManufacturing\\HybridManufacturing\\models\\coral\\cont_red_number.txt");
	ofstream ofile_cont_size("E:\\Hybrid manufacturing\\HybridManufacturing\\HybridManufacturing\\models\\coral\\cont_size.txt");
	/*bool* flag_sample_point_used = new bool[sampling.sample_points.size()];
	for (int i = 0; i < sampling.sample_points.size(); i++)
		flag_sample_point_used[i] = false;*/
	cout << "--------------------- height of outer-beam search: 1 ---------------------";
	while (last_step_nodes.size() != 0 || candidate_nodes.size() != 0) {
		/*if (height_of_beam_search ==2)
			W1 = 5;
		else
			W1 = 1;*/
		vector<vector<bool>> is_fragile_V_2;
		int now_last_node = -1;
		const bool current_node_loaded = PrepareOuterBeamSearchNode(
			last_step_nodes,
			is_fragile_V_2,
			now_last_node,
			height_of_beam_search,
			cont_number_of_queue,
			tree_entries,
			sum_time_5);
		//////////////

		if (current_node_loaded) {
			vector<bool> judge_S_be_searched;
			vector<bool> judge_covering_points_be_searched;
			sum_time_3 += UpdateSubtractiveDependencyGraph(
				ori_all_the_area_S,
				is_fragile_V_2,
				now_last_node,
				height_of_beam_search,
				cont_number_of_queue,
				tree_entries,
				judge_S_be_searched,
				judge_covering_points_be_searched);
			/////////////////////////////////////

			const double accessibility_marker_radius = std::max(0.25, dh * 0.15);
			const std::string accessibility_node_tag =
				MakeAccessibilityDebugNodeTag(height_of_beam_search, cont_number_of_queue, now_last_node);
			accessibility_visualization::WriteSubtractiveAccessibilityDebugVisualizations(
				UnaccessibleVisDir(),
				accessibility_node_tag,
				accessibility_marker_radius,
				V,
				judge_S_be_searched,
				judge_covering_points_be_searched,
				map_S_and_vertex,
				map_covering_points_and_vertex);
			if (!sampling_subtractive.sample_points.empty()) {
				const std::vector<bool> active_original_vertices =
					accessibility_visualization::BuildActiveOriginalVertexMask(
						V,
						Katana::Instance().vertices);
				int current_selected_s_id = -1;
				int current_selected_cell_id = -1;
				accessibility_visualization::WriteHighestZInaccessiblePointAllOrientationToolCollisionObj(
					UnaccessibleVisPath(
						"access_debug_subtractive_current_node_tool_collision"
						+ accessibility_node_tag
						+ ".obj"),
					V,
					all_voronoi_cells,
					sampling_subtractive.sample_points,
					map_S_and_vertex,
					cutting_tool,
					current_selected_s_id,
					current_selected_cell_id,
					&active_original_vertices,
					&judge_S_be_searched);
			}

			/*int end;
			if (height_of_beam_search != 4)
				end = sampling.sample_points.size();
			else
				end = 1;*/
			std::cout << sampling.sample_points.size() << endl;
			for (int ori = 0; ori < sampling.sample_points.size(); ori++) {	//枚举所有采样方向，作为增材分层方向
				printf("\r[%d%%]>", ori * 100 / (sampling.sample_points.size() - 1));
				for (int j = 1; j <= ori * 20 / sampling.sample_points.size(); j++)
					std::cout << "▉";

				Eigen::Vector3d vectorAfter;
				if (!PrepareOrientationSliceData(sampling, ori, now_last_node, save_ori, rotMatrix, vectorAfter)) {
					continue;
				}

				double temp_time_1 = 0;
				double temp_time_2 = 0;
				Layer_Graph layer_graph = BuildAdditiveLayerGraphWithSurfaceMesh(
					vectorAfter,
					height_of_beam_search,
					tree_entries[now_last_node].continue_id,
					the_nozzle,
					temp_time_1,
					temp_time_2);

				if (open_vis_additive_accessibility_debug) {
					accessibility_visualization::WriteAdditiveAccessibilityDebugVisualization(
						UnaccessibleVisDir(),
						layer_graph,
						accessibility_node_tag,
						ori,
						accessibility_marker_radius);
				}
				if (open_vis_additive_self_support_debug) {
					accessibility_visualization::WriteAdditiveRootLayerSelfSupportDebugVisualization(
						UnaccessibleVisDir(),
						layer_graph,
						accessibility_node_tag,
						ori,
						vectorAfter,
						current_node_mesh_rotated);
				}

				//Layer_Graph layer_graph = BuildAdditiveLayerGraph(
				//	vectorAfter,
				//	height_of_beam_search,
				//	tree_entries[now_last_node].continue_id,
				//	the_nozzle,
				//	sum_time_4,
				//	sum_time_2);

				//////////////////////////////////////
				vector<Eigen::MatrixXd> fragile_V;
				for (int i = 0; i < is_fragile_V_2[ori].size(); i++)
					if (is_fragile_V_2[ori][i] == true)
						fragile_V.push_back(V_2[i]);

				all_solutions_of_selected_layers.clear();
				all_solutions_of_selected_layers_contain.clear();
				all_solutions_of_selected_layer_holes.clear();
				all_cut_layers.clear(); all_cut_layers_dependency_layer.clear();
				pathes_include_S.clear(); pathes_include_sample_points.clear(); paths_include_covering_points.clear();
				//std::cout << "cc" << endl;


				bool flag_continue = false;	//是否继续从当前节点向下搜索
				bool jud_admit = true;
				if (!ProcessOrientationSearch(
					layer_graph,
					judge_S_be_searched,
					judge_covering_points_be_searched,
					W1,
					sum_time,
					flag_continue,
					jud_admit)) {
					continue;
				}

				final_pathes_include_S.insert(final_pathes_include_S.end(), pathes_include_S.begin(), pathes_include_S.end());
				final_pathes_include_sample_points.insert(final_pathes_include_sample_points.end(), pathes_include_sample_points.begin(), pathes_include_sample_points.end());
				final_pathes_include_covering_points.insert(final_pathes_include_covering_points.end(), paths_include_covering_points.begin(), paths_include_covering_points.end());

				for (int i = 0; i < all_solutions_of_selected_layers.size(); i++) {	//将当前方向生成的所有解加入树结构中，作为下一层搜索的节点
					index_node++;
					BeamTreeEntry entry;
					entry.layers = all_solutions_of_selected_layers[i];
					entry.contain = all_solutions_of_selected_layers_contain[i];
					if (i < static_cast<int>(all_solutions_of_selected_layer_holes.size())) {
						entry.holes = all_solutions_of_selected_layer_holes[i];
					}
					entry.cut_layers = all_cut_layers[i];
					entry.cut_layers_dependency_layer = all_cut_layers_dependency_layer[i];
					entry.fragile_v = fragile_V;
					entry.judge_s_be_searched = judge_S_be_searched;
					entry.parent_id = now_last_node;
					entry.pre_queue_index = cont_number_of_queue;
					candidate_nodes.push_back(index_node);
					vector<Eigen::Vector3d> temp_vec;
					temp_vec.push_back(vectorAfter);
					save_ori.push_back(temp_vec);
					entry.judge_continue = flag_continue;
					if (tree_entries[now_last_node].judge_continue == true)
						entry.continue_id = tree_entries[now_last_node].continue_id + 1;
					else
						entry.continue_id = 0;
					entry.error = false;
					tree_entries.push_back(entry);
				}
			}
		}



		if (current_node_loaded && tree_entries[now_last_node].judge_continue == false)
			last_step_nodes.pop();
		if (!current_node_loaded && !last_step_nodes.empty())
			last_step_nodes.pop();
		cont_number_of_queue++;
		if (last_step_nodes.size() == 0 || (current_node_loaded && tree_entries[now_last_node].judge_continue == true)) {
			cont_number_of_queue = 0;
			height_of_beam_search++;
			int cont_w = 0;  //应为0开始
			cout << endl << "Decomposing and sorting......" << endl;
			if (candidate_nodes.empty()) {
				std::cout << "[BeamSearch] No candidate nodes remain after skipping invalid mesh node(s). Stop outer beam search." << std::endl;
				break;
			}
			/////////sort_candidate_nodes//////////
			start_time = clock();
			vector<all_value> all_calculated_value;
			vector<all_value> pure_value;
			all_calculated_value.resize(candidate_nodes.size());
			pure_value.resize(candidate_nodes.size());
			for (int i = 0; i < static_cast<int>(candidate_nodes.size()); i++) {
				EvaluateCandidateNode(
					i,
					candidate_nodes,
					all_calculated_value,
					tree_entries,
					save_ori,
					height_of_beam_search,
					cont_number_of_queue,
					file_name,
					now_last_node,
					the_nozzle);
			}

			sort_candidate_nodes(candidate_nodes, tree_entries, final_pathes_include_S, all_calculated_value, final_pathes_include_covering_points, height_of_beam_search, save_ori, pure_value, tree_entries[candidate_nodes[0]].continue_id);
			end_time = clock();

			ofile_cont << final_pathes_include_S[candidate_nodes[0]].size() << endl;
			double temp_record_size = 0;
			for (int j = 0; j < tree_entries[candidate_nodes[0]].layers.size(); j++) {
				for (int k = 0; k < tree_entries[candidate_nodes[0]].layers[j].size() - 1; k++) {
					temp_record_size += distance3d(tree_entries[candidate_nodes[0]].layers[j][k], tree_entries[candidate_nodes[0]].layers[j][k + 1]);
				}
			}
			ofile_cont_size << temp_record_size << endl;
			//////////////////////////////////////
			if (tree_entries[now_last_node].judge_continue == true) {
				W1 = 1;
			}

			//cout << "bbb";

			bool jud_continue_last_node = false;
			vector<int> current_selected_nodes;
			//	for (int i = 0; i < candidate_nodes.size(); i++)
			//		cout << save_ori[candidate_nodes[i]][0].x() << " " << save_ori[candidate_nodes[i]][0].y() << " " << save_ori[candidate_nodes[i]][0].z() << " " << endl;
			start_time_7 = clock();
			int accepted_candidate_count = 0;
			const int candidate_attempt_limit = tree_entries[now_last_node].judge_continue
				? std::min(W1, static_cast<int>(candidate_nodes.size()))
				: static_cast<int>(candidate_nodes.size());
			while (cont_w < candidate_attempt_limit && accepted_candidate_count < W1) {
				const int selected_node = candidate_nodes[cont_w];
				//decompose the model for every node//
				int index_of_pre_node = tree_entries[selected_node].pre_queue_index;
				tree_entries[selected_node].source_input_file = make_cut_source_file(
					tree_entries[now_last_node].judge_continue,
					height_of_beam_search,
					index_of_pre_node,
					tree_entries[selected_node].continue_id);
				vector<vector<cv::Point3d>> all_cut_layers;
				vector<int> flag_cut_layers_is_hole;
				vector<int> all_cut_layers_dependency_layer;
				BuildCutLayersForCandidate(
					cont_w,
					candidate_nodes,
					tree_entries,
					all_cut_layers,
					all_cut_layers_dependency_layer,
					flag_cut_layers_is_hole);

				cout << "***** number_of_S: " << final_pathes_include_S[candidate_nodes[cont_w]].size() << endl;
				cout << "***** all_cut_layers.size(): " << all_cut_layers.size() << endl;

				std::vector<TRiangle> current_remove_triangles;
				Slicer_2 current_slicer;
				//cout << "aaa";
				/*for (int i = 0; i < Tree_nodes_continue_id.size(); i++)
					cout << "*"<<Tree_nodes_continue_id[i] << endl;*/

				CutMesh(
					tree_entries[selected_node].layers,
					tree_entries[selected_node].contain,
					tree_entries[selected_node].holes,
					all_cut_layers, save_ori[selected_node][0],
					height_of_beam_search,
					cont_number_of_queue,
					index_of_pre_node,
					all_cut_layers_dependency_layer,
					jud_outer_beam_search_terminate,
					current_remove_triangles,
					current_slicer,
					tree_entries[selected_node].judge_continue, tree_entries[now_last_node].judge_continue,
					tree_entries[now_last_node].pre_queue_index,
					tree_entries[selected_node].error,
					selected_node,
					tree_entries[selected_node].continue_id,
					flag_cut_layers_is_hole,
					tree_entries[selected_node].judge_s_be_searched);

				if (tree_entries[selected_node].error) {
					std::cerr << "[BeamSearch] Skip invalid cut candidate: "
						<< "node=" << selected_node
						<< ", height=" << height_of_beam_search
						<< ", queue_index=" << cont_number_of_queue
						<< std::endl;
					++cont_w;
					continue;
				}

				LogSelectedCandidateMetrics(
					cont_w,
					candidate_nodes,
					now_last_node,
					save_ori,
					tree_entries,
					last_step_nodes,
					all_calculated_value,
					pure_value,
					Sum_candidate_blocks,
					Sum_connected_components,
					cont_extra_additive_orientation,
					evaluation_value);


				//if (Tree_nodes_judge_continue[now_last_node] == true) {
				//	cout << "aaa";
				//	Tree_nodes.pop_back();
				//	Tree_nodes_cut_layers.pop_back();
				//	Tree_nodes_num_of_cut_layers_dependency_layer.pop_back();
				//	Tree_nodes_fragile_V.pop_back();
				//	index_node--;
				//	candidate_nodes.pop_back();
				//	save_ori.pop_back();
				//	Tree_nodes_judge_continue.pop_back();
				//	Tree_nodes_error.pop_back();
				//	height_of_beam_search--;
				//	Tree_nodes_judge_continue[now_last_node] = false;
				//	//Tree_nodes_continue_id[now_last_node]++;
				//}
				cout << "&&&&" << height_of_beam_search << endl;
				if (tree_entries[selected_node].judge_continue == true) {
					cout << "AAAA!" << endl;
					jud_continue_last_node = true;

				}
				if (tree_entries[now_last_node].judge_continue == true) { //Tree_nodes_judge_continue[now_last_node] == true  //应该改为Tree_nodes_judge_continue[candidate_nodes[0]] == true?
					cout << "BBBB!" << candidate_nodes[0] << endl;
					if (tree_entries[candidate_nodes[0]].judge_continue == false)
						tree_entries[now_last_node].judge_continue = false;
					height_of_beam_search--;
					tree_entries.pop_back();
					index_node--;
					candidate_nodes.pop_back();
					//Tree_nodes_judge_continue[now_last_node] = false;
					tree_entries[now_last_node].continue_id = tree_entries[candidate_nodes[0]].continue_id;
					save_ori[now_last_node].push_back(save_ori[candidate_nodes[0]][0]);
					save_ori.pop_back();

				}
				if (selected_node >= 0 && selected_node < static_cast<int>(tree_entries.size())) {
					current_selected_nodes.push_back(selected_node);
				}


				int subtractive_accessibility_patch_count = 0;
				if (height_of_beam_search <= 150) {
					//cout << save_ori[candidate_nodes[cont_w]].x() <<" "<< save_ori[candidate_nodes[cont_w]].y() << " " << save_ori[candidate_nodes[cont_w]].z() << endl;
					subtractive_accessibility_patch_count = subtractive_accessibility_decomposition(current_remove_triangles, height_of_beam_search, cont_number_of_queue, cutting_tool, current_slicer);
					//subtractive_accessibility_decomposition(current_remove_triangles, height_of_beam_search, index_of_pre_node, cutting_tool, current_slicer);
				}

				subtractive_remove_output(current_remove_triangles, current_slicer, height_of_beam_search, cont_number_of_queue);
				if (selected_node >= 0 && selected_node < static_cast<int>(tree_entries.size())) {
					tree_entries[selected_node].subtractive_accessibility_patch_count = subtractive_accessibility_patch_count;
					tree_entries[selected_node].patch_output_file =
						PatchVisPath(
							"block_patch-" + to_string(height_of_beam_search)
							+ "_" + to_string(cont_number_of_queue)
							+ ".obj");
				}

				cont_number_of_queue++;
				//////////////////////////////////////

				accepted_candidate_count++;
				cont_w++;
			}
			if (!current_selected_nodes.empty()) {
				last_selected_nodes = current_selected_nodes;
			}
			end_time_7 = clock();
			std::cout << "&&&time&&& Evaluation " << double(end_time - start_time + end_time_7 - start_time_7) / CLOCKS_PER_SEC << std::endl;
			std::cout << "&&&time&&& Update subtractive dependency graph: " << sum_time_3 << std::endl;
			std::cout << "&&&time&&& Slicing: " << sum_time_4 << std::endl;
			std::cout << "&&&time&&& Build addictive dependency graph: " << sum_time_2 << std::endl;
			std::cout << "&&&time&&& Co graph merging: " << sum_time << std::endl;
			total_time_1 += sum_time_3;
			total_time_2 += sum_time_2 + sum_time_4 + sum_time_5;
			total_time_3 += sum_time;
			total_time_4 += double(end_time - start_time + end_time_7 - start_time_7) / CLOCKS_PER_SEC;
			sum_time = 0;
			sum_time_2 = 0;
			sum_time_3 = 0;
			sum_time_4 = 0;
			sum_time_5 = 0;
			if (accepted_candidate_count == 0) {
				std::cout << "[BeamSearch] All attempted candidates at height "
					<< height_of_beam_search
					<< " were rejected; no next-layer node was generated. Stop outer beam search."
					<< std::endl;
				candidate_nodes.clear();
				break;
			}
			//std::cout << "&&&time&&& a step of beam search: " << double(end_time_2 - start_time_2) / CLOCKS_PER_SEC << std::endl;
			/*if (height_of_beam_search == 5)
				W1 = 20;*/
				/*if (height_of_beam_search >= 2 && Tree_nodes_continue_id[candidate_nodes[0]] >= 0) {
					end_time_total = clock();
					std::cout << "&&&time&&& total of outer-beam search: " << double(end_time_total - start_time_total) / CLOCKS_PER_SEC << "s &&&&&&&" << std::endl << endl;
					system("pause");
				} */

				//system("pause");
				//save_ori.clear();
			cont_number_of_queue = 0;
			if (jud_outer_beam_search_terminate == true) {
				/*cout << "$$$$$$$$$$$$$$$$$$$$$$$$$ Data $$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
				cout << "#S: " << all_voronoi_cells.size() << endl;
				cout << "#P_i: " << num_inaccessible_points << endl;
				cout << "#R_p: " << float(num_inaccessible_points) / float(all_voronoi_cells.size()) << endl;
				cout << "#B_C !!: " << Sum_candidate_blocks << endl;
				cout << "#I: " << height_of_beam_search << endl;
				cout << "#B: " << height_of_beam_search << endl;
				cout << "#O: " << cont_extra_additive_orientation + height_of_beam_search << endl;
				cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;

				cout << "$$$$$$$$$$$$$$$$$$$$$$$$$ Total time $$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
				cout << "Build subtractive graph: " << time_build_subtractive_graph << endl;
				cout << "Update subtractive dependency graph: " << total_time_1 << endl;
				cout << "Build addictive dependency graph: " << total_time_2 << endl;
				cout << "Co graph merging: " << total_time_3 << endl;
				cout << "Evaluation: " << total_time_4 << endl;
				cout << "Total: " << time_build_subtractive_graph + total_time_1 + total_time_2 + total_time_3 + total_time_4 << endl;
				cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;*/

				double global_score = 1 / (height_of_beam_search);
				fstream f;
				//追加写入,在原来基础上加了ios::app  faith
				std::string txt_name = file_name + ".txt";
				f.open(txt_name.c_str(), ios::out | ios::app);
				//输入你想写入的内容 
				f << global_score << endl;
				f.close();

				break;
			}
			if (height_of_beam_search >= 10) {
				fstream f;
				//追加写入,在原来基础上加了ios::app  faith
				std::string txt_name = file_name + ".txt";
				f.open(txt_name.c_str(), ios::out | ios::app);
				//输入你想写入的内容 
				f << 0 << endl;
				f.close();
				//break;
			}
			if (jud_continue_last_node == true)
				cout << endl << "--------------------- height of outer-beam search: " << height_of_beam_search << " ---------------------";
			else
				cout << endl << "--------------------- height of outer-beam search: " << height_of_beam_search + 1 << " ---------------------";
			candidate_nodes.clear();
		}
	}

	const std::string ancestor_source_report_file =
		AncestorSourceVisPath(
			"beam_search_last_layer_ancestor_sources_"
			+ MakeAccessibilityDebugFileToken(file_name)
			+ ".txt");
	std::ofstream ancestor_source_report(ancestor_source_report_file);
	std::ostream& ancestor_out = ancestor_source_report.is_open() ? ancestor_source_report : std::cout;
	if (!ancestor_source_report.is_open()) {
		std::cout << "[BeamSearch] Failed to open ancestor source report: "
			<< ancestor_source_report_file << std::endl;
	}

	ancestor_out << "[BeamSearch] Last layer ancestor source files\n";
	ancestor_out << std::setprecision(17);
	ancestor_out << "last_selected_node_count: " << last_selected_nodes.size() << '\n';
	for (int final_node : last_selected_nodes) {
		if (final_node < 0 || final_node >= static_cast<int>(tree_entries.size())) {
			continue;
		}

		ancestor_out << "final_node: " << final_node << '\n';
		int node = final_node;
		int depth = 0;
		while (node >= 0 && node < static_cast<int>(tree_entries.size())) {
			const auto& entry = tree_entries[node];
			const std::string source_file = entry.source_input_file.empty()
				? (file_name + "-0_0.obj")
				: entry.source_input_file;
			const bool has_saved_orientations =
				node > 0
				&& node < static_cast<int>(save_ori.size())
				&& !save_ori[node].empty();
			ancestor_out << "  depth_from_final: " << depth
				<< ", node: " << node
				<< ", parent: " << entry.parent_id
				<< ", pre_queue_index: " << entry.pre_queue_index
				<< ", continue_id: " << entry.continue_id
				<< ", source_input_file: " << source_file
				<< ", subtractive_accessibility_patch_count: " << entry.subtractive_accessibility_patch_count
				<< ", patch_output_file: " << (entry.patch_output_file.empty() ? "<none>" : entry.patch_output_file)
				<< ", additive_orientation_count: "
				<< (has_saved_orientations ? save_ori[node].size() : 0)
				<< ", additive_orientations: ";

			if (node == 0) {
				ancestor_out << "<root>";
			}
			else if (!has_saved_orientations) {
				ancestor_out << "<missing>";
			}
			else {
				ancestor_out << '[';
				for (std::size_t orientation_index = 0;
					orientation_index < save_ori[node].size();
					++orientation_index) {
					if (orientation_index != 0) {
						ancestor_out << "; ";
					}
					const Eigen::Vector3d& orientation = save_ori[node][orientation_index];
					ancestor_out << '(' << orientation.x()
						<< ", " << orientation.y()
						<< ", " << orientation.z() << ')';
				}
				ancestor_out << ']';
			}
			ancestor_out << '\n';

			if (entry.parent_id == node) {
				ancestor_out << "  stop: parent_id points to itself\n";
				break;
			}
			node = entry.parent_id;
			++depth;
		}
	}
	if (ancestor_source_report.is_open()) {
		ancestor_source_report.close();
		std::cout << "[BeamSearch] Last layer ancestor source report saved to "
			<< ancestor_source_report_file << std::endl;
	}

	vector<int> last_available_block;
	vector<int> exist_points;
	vector<vector<int>> is_point_exist_in_block;
	vector<vector<int>> is_point_available_in_block;
	vector<int> is_available(height_of_beam_search + 1, 1);


	last_available_block.resize(V.rows());
	for (int i = 0; i < V.rows(); i++) {
		last_available_block[i] = height_of_beam_search;
	}

	for (int i = 0; i < V.rows(); i++) {
		is_point_exist_in_block.push_back(is_available);
	}

	for (int i = 0; i < V.rows(); i++)
		is_point_available_in_block.push_back(is_available);

	cout << "TTYTYTYTYTY";

	//subtractive_accessibility_decomposition_global(height_of_beam_search, cutting_tool);
}

double HybridManufacturing::UpdateSubtractiveDependencyGraph(
	const vector<vector<area_S>>& ori_all_the_area_S,
	vector<vector<bool>>& is_fragile_V_2,
	int now_last_node,
	int height_of_beam_search,
	int cont_number_of_queue,
	const vector<BeamTreeEntry>& tree_entries,
	vector<bool>& judge_S_be_searched,
	vector<bool>& judge_covering_points_be_searched)
{
	clock_t start_time_5 = clock();
	vector<int> S_in_block;
	vector<int> covering_points_in_block;	//记录被碰撞单元索引
	judge_S_be_searched.clear();
	judge_covering_points_be_searched.clear();
	judge_S_be_searched.resize(all_the_area_S.size());
	judge_covering_points_be_searched.resize(all_the_covering_points.size());
	for (int i = 0; i < all_the_area_S.size(); i++)
		judge_S_be_searched[i] = false;
	for (int i = 0; i < all_the_covering_points.size(); i++)
		judge_covering_points_be_searched[i] = false;
	S_in_block.clear();
	covering_points_in_block.clear();
	for (int i = 0; i < V.rows(); i++) {	//枚举原始输入模型网格V的所有顶点，判断该顶点是否仍存在于当前块中
		bool jud_still_exist = false;
		/*if (flag_sample_point_used[i])
			continue;*/
		for (int j = 0; j < Katana::Instance().vertices.size(); j++) {
			if (abs(V(i, 0) - Katana::Instance().vertices[j].x) <= 0.001 && abs(V(i, 1) - Katana::Instance().vertices[j].y) <= 0.001 && abs(V(i, 2) - Katana::Instance().vertices[j].z) <= 0.001) {
				jud_still_exist = true;	//标记该顶点仍存在
				//flag_sample_point_used[i] = true;
				break;
			}
		}

		if (jud_still_exist == false) {	//如果该顶点不存在于当前块中
			for (int k = 0; k < is_fragile_V_2.size(); k++)
				is_fragile_V_2[k][i] = false;	//is_fragile_V_2中对应该顶点的信息全部置为false
			const auto area_s_it = map_S_and_vertex_inv.find(i);
			if (area_s_it != map_S_and_vertex_inv.end()) {
				S_in_block.push_back(area_s_it->second);
			}
			const auto covering_point_it = map_covering_points_and_vertex_inv.find(i);
			if (covering_point_it != map_covering_points_and_vertex_inv.end()) {
				covering_points_in_block.push_back(covering_point_it->second);
			}
		}
	}

	/*while (pre_tree_nodes[now_last_node] != -1) {
		for (int i = 0; i < final_pathes_include_S[now_last_node].size(); i++) {
			S_in_block.push_back(final_pathes_include_S[now_last_node][i]);
		}
		for (int i = 0; i < final_pathes_include_sample_points[now_last_node].size(); i++) {
			sample_points_in_block.push_back(final_pathes_include_sample_points[now_last_node][i]);
		}
		now_last_node = pre_tree_nodes[now_last_node];
	}
	now_last_node = last_step_nodes.front();*/

	for (int i = 0; i < S_in_block.size(); i++)
		judge_S_be_searched[S_in_block[i]] = true;

	ori_num_points_of_ori_in_all_the_area_S.assign(
		ori_all_the_area_S.size(),
		vector<int>(sampling_subtractive.sample_points.size(), 0)); //为每个不可达点在ori_num_points_of_ori_in_all_the_area_S分配一个向量
	//每个不可达点记录在各个采样方向上对应的点数

	for (size_t point_id = 0; point_id < ori_all_the_area_S.size(); ++point_id) {
		const auto& orientations_for_point = ori_all_the_area_S[point_id];
		auto& orientation_point_count = ori_num_points_of_ori_in_all_the_area_S[point_id];
		for (const auto& area : orientations_for_point) {
			++orientation_point_count[area.oriId];
		}
		for (int ori = 0; ori < sampling_subtractive.sample_points.size(); ++ori) {
			if (orientation_point_count[ori] == 0) {
				if (height_of_beam_search <= 1) {
					std::cout << "ori_num_points_of_ori_in_all_the_area_S[point_id][ori] is zero for point_id " << point_id << " and ori " << ori << std::endl;
				}
			}
		}
	}

	for (int covering_point_id : covering_points_in_block) {
		if (covering_point_id < 0
			|| covering_point_id >= static_cast<int>(all_the_covering_points.size())) {
			std::cout << "[UpdateSubtractiveDependencyGraph] Invalid covering point id: "
				<< covering_point_id << std::endl;
			continue;
		}
		const auto& covering_points_for_vertex = all_the_covering_points[covering_point_id];
		for (const auto& covering_point : covering_points_for_vertex) {
			if (covering_point.pointId < 0
				|| covering_point.pointId >= static_cast<int>(ori_num_points_of_ori_in_all_the_area_S.size())
				|| covering_point.oriId < 0
				|| covering_point.oriId >= static_cast<int>(ori_num_points_of_ori_in_all_the_area_S[covering_point.pointId].size())) {
				std::cout << "[UpdateSubtractiveDependencyGraph] Invalid covering relation: pointId="
					<< covering_point.pointId
					<< ", oriId=" << covering_point.oriId
					<< std::endl;
				continue;
			}
			ori_num_points_of_ori_in_all_the_area_S[covering_point.pointId][covering_point.oriId]--;
			int cover_count = ori_num_points_of_ori_in_all_the_area_S[covering_point.pointId][covering_point.oriId];
			if (cover_count < 0) {
				std::cout << "Warning: Negative count for pointId " << covering_point.pointId
					<< ", oriId " << covering_point.oriId
					<< ". Current count: " << cover_count << std::endl;
			}
			if (cover_count == 0) {
				judge_S_be_searched[covering_point.pointId] = true;
			}
		}
	}

	for (int covering_point_id : covering_points_in_block) {
		judge_covering_points_be_searched[covering_point_id] = true;
	}

	if (tree_entries[now_last_node].judge_continue == true) {
		Katana::Instance().stl.loadStl((file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_" + to_string(tree_entries[now_last_node].continue_id) + "_subblock.stl").c_str());
		Katana::Instance().temp_vertices.resize(Katana::Instance().vertices.size());
		Katana::Instance().temp_vertices = Katana::Instance().vertices;
	}
	clock_t end_time_5 = clock();
	return double(end_time_5 - start_time_5) / CLOCKS_PER_SEC;
}

bool HybridManufacturing::PrepareOrientationSliceData(
	const SAMPLE_ON_BALL& sampling,
	int ori,
	int now_last_node,
	const vector<vector<Eigen::Vector3d>>& save_ori,
	Eigen::Matrix3d& rot_matrix,
	Eigen::Vector3d& vector_after)
{
	Eigen::Vector3d vector_continue(sampling.sample_points[ori]);  //与祖先Node的ori相同时容易出现问题，暂时先避免用相同ori
	bool jud_continue = false;

	for (int j = 0; j < save_ori[now_last_node].size(); j++)
		if (vector_continue == save_ori[now_last_node][j]) {
			jud_continue = true;
			break;
		}

	if (jud_continue == true && ori != 0) {
		return false;
	}

	std::vector<Eigen::MatrixXd> temp_V;	//记录当前块的顶点信息，用于旋转
	temp_V.resize(Katana::Instance().temp_vertices.size());

	V_2.resize(V.rows());	//V_2记录输入模型网格V的顶点信息，用于旋转
	for (int i = 0; i < Katana::Instance().temp_vertices.size(); i++) {
		temp_V[i].resize(3, 1);
		temp_V[i](0, 0) = Katana::Instance().temp_vertices[i].x;
		temp_V[i](1, 0) = Katana::Instance().temp_vertices[i].y;
		temp_V[i](2, 0) = Katana::Instance().temp_vertices[i].z;
	}
	for (int i = 0; i < V.rows(); i++) {
		V_2[i].resize(3, 1);
		V_2[i](0, 0) = V.row(i).x();
		V_2[i](1, 0) = V.row(i).y();
		V_2[i](2, 0) = V.row(i).z();
	}
	Eigen::Vector3d vector_before(0, 0, 1);
	vector_after = Eigen::Vector3d(sampling.sample_points[ori]);
	rot_matrix = Eigen::Quaterniond::FromTwoVectors(vector_before, vector_after).toRotationMatrix();
	for (int i = 0; i < Katana::Instance().vertices.size(); i++)
		temp_V[i] = rot_matrix.inverse() * temp_V[i];	//将temp_V里的顶点坐标旋转，使得增材打印方向与z轴平行
	for (int i = 0; i < V.rows(); i++)
		V_2[i] = rot_matrix.inverse() * V_2[i];	//将V_2里的顶点坐标旋转，使得增材打印方向与z轴平行
	for (int i = 0; i < Katana::Instance().vertices.size(); i++) {
		Katana::Instance().vertices[i].x = temp_V[i](0, 0);
		Katana::Instance().vertices[i].y = temp_V[i](1, 0);
		Katana::Instance().vertices[i].z = temp_V[i](2, 0);	//将temp_V里的顶点坐标更新回Katana的vertices
	}

	for (int i = 0; i < Katana::Instance().triangles.size(); i++) {	//将每个三角形的顶点按z值从小到大排序，便于后续分层处理
		if (Katana::Instance().triangles[i].vertices[0]->z > Katana::Instance().triangles[i].vertices[1]->z) std::swap(Katana::Instance().triangles[i].vertices[0], Katana::Instance().triangles[i].vertices[1]);
		if (Katana::Instance().triangles[i].vertices[0]->z > Katana::Instance().triangles[i].vertices[2]->z) std::swap(Katana::Instance().triangles[i].vertices[0], Katana::Instance().triangles[i].vertices[2]);
		if (Katana::Instance().triangles[i].vertices[1]->z > Katana::Instance().triangles[i].vertices[2]->z) std::swap(Katana::Instance().triangles[i].vertices[1], Katana::Instance().triangles[i].vertices[2]);
	}

	current_node_mesh_rotated.clear();
	current_node_mesh_rotated = current_node_mesh;

	for (auto it : current_node_mesh_rotated.vertices())
	{
		auto& vertex = current_node_mesh_rotated.point(it);
		auto temp_vec = Eigen::Vector3d(vertex.x(), vertex.y(), vertex.z());
		temp_vec = rot_matrix.inverse() * temp_vec;
		vertex = Point_3(temp_vec.x(), temp_vec.y(), temp_vec.z());
	}

	CGAL::IO::write_polygon_mesh("model\\temp_current_node_mesh.obj", current_node_mesh);
	CGAL::IO::write_polygon_mesh("model\\temp_rotated.obj", current_node_mesh_rotated);

	return true;
}

// 保持函数签名与声明一致，避免不必要的声明差异。
Layer_Graph HybridManufacturing::BuildAdditiveLayerGraph(
	const Eigen::Vector3d& vector_after,
	int height_of_beam_search,
	int continue_node_id,
	const nozzle& the_nozzle,
	double& slicing_time,
	double& graph_time)
{
	clock_t start_time_6 = clock();
	vector<vector<vector<Vertex>>> all_slice_points; //(i,j,k) 第i层第j个loop的第k个点的坐标
	vector<vector<vector<Vertex>>> all_slice_points_contain;
	Katana::Instance().slicer.buildLayers();


	bool segments_manifold = Katana::Instance().slicer.buildSegments();
	if (!segments_manifold) {
		std::cout << "[HybridManufacturing::BuildAdditiveLayerGraph] !segments_manifold" << std::endl;
		//Visualize_layer_segments(Katana::Instance().layers);

		CGAL::Polygon_mesh_slicer<SurfaceMesh, Kernel> slicer_cgal(current_node_mesh_rotated);
		std::vector<Polylines> layer_polylines;
		layer_polylines.resize(Katana::Instance().layers.size());
		for (auto& lay : Katana::Instance().layers) {
			slicer_cgal(Kernel::Plane_3(0, 0, -1, lay.z), std::back_inserter(layer_polylines[&lay - &Katana::Instance().layers[0]]));
		}

		//Visualize_layer_polylines(layer_polylines);
	}

	Katana::Instance().gcode.write(all_slice_points, all_slice_points_contain);	//katana分层并得到边界轮廓

	vector<vector<vector<Eigen::Vector3d>>> vector_after_all_slice_points;	//记录每个切片轮廓上每个点的坐标
	vector<vector<vector<Eigen::Vector3d>>> vector_after_all_slice_points_contain;
	vector_after_all_slice_points.resize(all_slice_points.size());
	vector_after_all_slice_points_contain.resize(all_slice_points_contain.size());
	for (int i = 0; i < all_slice_points.size(); i++) {
		vector_after_all_slice_points[i].resize(all_slice_points[i].size());
		for (int j = 0; j < all_slice_points[i].size(); j++) {
			vector_after_all_slice_points[i][j].resize(all_slice_points[i][j].size());
			for (int k = 0; k < all_slice_points[i][j].size(); k++) {
				vector_after_all_slice_points[i][j][k] = Eigen::Vector3d(all_slice_points[i][j][k].x, all_slice_points[i][j][k].y, all_slice_points[i][j][k].z);
			}
		}
	}
	for (int i = 0; i < all_slice_points_contain.size(); i++) {
		vector_after_all_slice_points_contain[i].resize(all_slice_points_contain[i].size());
		for (int j = 0; j < all_slice_points_contain[i].size(); j++) {
			vector_after_all_slice_points_contain[i][j].resize(all_slice_points_contain[i][j].size());
			for (int k = 0; k < all_slice_points_contain[i][j].size(); k++) {
				vector_after_all_slice_points_contain[i][j][k] = Eigen::Vector3d(all_slice_points_contain[i][j][k].x, all_slice_points_contain[i][j][k].y, all_slice_points_contain[i][j][k].z);
			}
		}
	}

	clock_t end_time_6 = clock();
	slicing_time += double(end_time_6 - start_time_6) / CLOCKS_PER_SEC;
	/////////////////////////////////////
	//std::cout << "aa" << Katana::Instance().vertices.size()<<" " << all_slice_points.size() << endl;
	//generate additive dependency graph//
	std::vector<Data> data;
	data.resize(1);
	data[0].ReadData(vector_after_all_slice_points, vector_after_all_slice_points_contain);
	Layer_Graph layer_graph(data[0]);
	clock_t start_time_4 = clock();
	layer_graph.GetTrianglesForLayers(all_slice_points, Katana::Instance().map_segment_triangles, Katana::Instance().vertices, vector_after, height_of_beam_search, continue_node_id);	//将切片轮廓映射到三角形集，建立每层三角形集合等中间信息？
	layer_graph.GenerateDependencyEdges();	//生成增材分层依赖图的边
	layer_graph.CollisionDetectionForAdditiveManufacturing(the_nozzle);	//增材的碰撞检测，标记增材的不可达点等信息
	clock_t end_time_4 = clock();
	graph_time += double(end_time_4 - start_time_4) / CLOCKS_PER_SEC;

	return layer_graph;
}

bool HybridManufacturing::ProcessOrientationSearch(
	Layer_Graph& layer_graph,
	const vector<bool>& judge_S_be_searched,
	const vector<bool>& judge_covering_points_be_searched,
	int& W1,
	double& sum_time,
	bool& flag_continue,
	bool& jud_admit)
{
	bool previous_is_continue = false;	//记录上一个节点是否为continue节点
	clock_t start_time_3 = clock();

	if (open_change_orientation == true)
		W1 = 1;
	DFS_search(layer_graph, flag_continue, previous_is_continue, judge_S_be_searched, judge_covering_points_be_searched, jud_admit);	//对当前方向进行增材分层依赖图的深度优先搜索，生成这个分层方向的解

	/*if (height_of_beam_search==1&& Tree_nodes_continue_id[now_last_node] >= 2)
		flag_continue = false;*/

	clock_t end_time_3 = clock();
	if (jud_admit == false)
		return false;
	sum_time += double(end_time_3 - start_time_3) / CLOCKS_PER_SEC;

	//cout << "()()()(" << double(end_time_3 - start_time_3) / CLOCKS_PER_SEC << endl;
	return true;
}

void HybridManufacturing::DFS_search(Layer_Graph layer_graph, bool& flag_continue, bool previous_is_continue, vector<bool> judge_S_be_searched, vector<bool> judge_covering_points_be_searched, bool& jud_admit)
{
	int W2 = 1;  //2
	struct DfsTreeEntry {
		int layer = -1;
		vector<int> area_s;
		vector<int> sample_points;
		vector<int> covering_points;
		int parent = -1;
	};
	vector<DfsTreeEntry> tree_entries;
	vector<int> candidate_nodes;
	queue<int> last_step_nodes;
	vector<int> S_in_block;
	vector<int> terminate_nodes;
	DfsTreeEntry root_entry;
	tree_entries.push_back(root_entry);
	int index_node = 0;
	last_step_nodes.push(index_node);
	//vector<vector<area_S>> ori_all_the_area_S = all_the_area_S;

	//std::cout << "AA" << endl;
	//find boundary point of every layer//
	vector<double> left_point_of_layer, right_point_of_layer, top_point_of_layer, bottom_point_of_layer;
	left_point_of_layer.resize(layer_graph.total_node_num); right_point_of_layer.resize(layer_graph.total_node_num);
	top_point_of_layer.resize(layer_graph.total_node_num); bottom_point_of_layer.resize(layer_graph.total_node_num);

	const double layer_polygon_offset = std::max(1.0, dh * 0.02);
	const double layer_polygon_repair_eps = std::max(1e-9, layer_polygon_offset * 0.01);
	const double layer_z_slack = std::max(1e-4, dh * 0.05);
	vector<vector<Point_2>> vec_points_2d(layer_graph.total_node_num);
	vector<LayerContainmentPolygon> vec_polygon(layer_graph.total_node_num);
	vector<vector<LayerContainmentPolygon>> vec_hole_polygons(layer_graph.total_node_num);
	vector<double> vec_polygon_z(layer_graph.total_node_num, 0.0);
	/*Point_2 pp(V(i, 0), V(i, 1));
	if (polygon.bounded_side(pp) != CGAL::ON_UNBOUNDED_SIDE)
		jud_still_exist = true;*/

	for (int i = 0; i < layer_graph.total_node_num; i++) {
		left_point_of_layer[i] = MAX_I;
		right_point_of_layer[i] = -MAX_I;
		top_point_of_layer[i] = -MAX_I;
		bottom_point_of_layer[i] = MAX_I;
	}
	int temp_num = 0;
	for (int i = 0; i < layer_graph.data.slice_points.size(); i++) {
		for (int j = 0; j < layer_graph.data.slice_points[i].size(); j++) {
			for (int k = 0; k < layer_graph.data.slice_points[i][j].size(); k++) {
				/*if (layer_graph.data.slice_points[i][j][k].x < left_point_of_layer[temp_num])
					left_point_of_layer[temp_num] = layer_graph.data.slice_points[i][j][k].x;
				if (layer_graph.data.slice_points[i][j][k].x > right_point_of_layer[temp_num])
					right_point_of_layer[temp_num] = layer_graph.data.slice_points[i][j][k].x;*/
					/*if (layer_graph.data.slice_points[i][j][k].y < bottom_point_of_layer[temp_num])
						bottom_point_of_layer[temp_num] = layer_graph.data.slice_points[i][j][k].y;
					if (layer_graph.data.slice_points[i][j][k].y > top_point_of_layer[temp_num])
						top_point_of_layer[temp_num] = layer_graph.data.slice_points[i][j][k].y;*/
				vec_points_2d[temp_num].push_back(Point_2(layer_graph.data.slice_points[i][j][k].x(), layer_graph.data.slice_points[i][j][k].y()));
			}
			if (!layer_graph.data.z_value[i][j].empty()) {
				vec_polygon_z[temp_num] = layer_graph.data.z_value[i][j][0];
			}
			vec_polygon[temp_num] = PrepareLayerContainmentPolygon(vec_points_2d[temp_num], layer_polygon_repair_eps);
			if (i < static_cast<int>(layer_graph.data.slice_points_holes.size())
				&& j < static_cast<int>(layer_graph.data.slice_points_holes[i].size())) {
				for (const auto& hole : layer_graph.data.slice_points_holes[i][j]) {
					vector<Point_2> hole_points;
					hole_points.reserve(hole.size());
					for (const Eigen::Vector2d& point : hole) {
						hole_points.emplace_back(point.x(), point.y());
					}
					vec_hole_polygons[temp_num].push_back(
						PrepareLayerContainmentPolygon(hole_points, layer_polygon_repair_eps));
				}
			}
			if (!vec_polygon[temp_num].is_simple) {
				std::cout << "[HybridManufacturing::DFS_search] layer polygon is not simple after local repair: layer="
					<< temp_num
					<< ", raw_points=" << vec_points_2d[temp_num].size()
					<< ", repaired_points=" << vec_polygon[temp_num].points.size()
					<< ", use even-odd containment with boundary offset=" << layer_polygon_offset
					<< std::endl;

				const std::string poly_file_name = VisPath("dfs_layer_not_simple_" + std::to_string(temp_num) + ".txt");
				std::ofstream poly_file(poly_file_name);
				if (poly_file.is_open()) {
					poly_file << std::setprecision(17);
					poly_file << "raw_point_count: " << vec_points_2d[temp_num].size() << '\n';
					poly_file << "repaired_point_count: " << vec_polygon[temp_num].points.size() << '\n';
					poly_file << "offset: " << layer_polygon_offset << '\n';
					poly_file << "raw_points:\n";
					for (const auto& point : vec_points_2d[temp_num]) {
						poly_file << CGAL::to_double(point.x()) << " " << CGAL::to_double(point.y()) << '\n';
					}
					poly_file << "repaired_points:\n";
					for (const auto& point : vec_polygon[temp_num].points) {
						poly_file << CGAL::to_double(point.x()) << " " << CGAL::to_double(point.y()) << '\n';
					}
				}
			}
			temp_num++;
		}
	}


	//////////////////////////////////////
	//std::cout << "BB" << endl;
	//find all the sample points in every layer//
	vector<vector<int>> sample_point_in_layer;	//存储每个layer中包含的sample point的索引
	sample_point_in_layer.resize(layer_graph.data.total_node_num);
	vector<bool> judge_sample_point_be_searched;
	judge_sample_point_be_searched.resize(V_2.size());
	for (int i = 0; i < V_2.size(); i++)
		judge_sample_point_be_searched[i] = false;
	for (int i = 0; i < layer_graph.data.total_node_num; i++) {
		for (int j = 0; j < V_2.size(); j++) {
			pair<int, int> index_slice_layer = layer_graph.data.index[i];
			Point_2 pp(V_2[j](0, 0), V_2[j](1, 0));
			if (judge_sample_point_be_searched[j] == false
				&& IsPointInsidePreparedMaterialRegion(
					vec_polygon[i], vec_hole_polygons[i], pp, layer_polygon_offset)
				&& IsPointWithinLayerZBand(
					layer_graph.data.z_value[index_slice_layer.first][index_slice_layer.second][0],
					V_2[j](2, 0),
					dh * 0.5,
					layer_z_slack)) {
				sample_point_in_layer[i].push_back(j);	//V_2的第j个点在第i个layer中
				judge_sample_point_be_searched[j] = true;	//标记该点需要被搜索
			}
		}
	}
	////////////////////////////////////////

	//std::cout << "CC" << endl;
	//find area S and covering points in the node//
	vector<vector<int>> temp_all_S_in_the_block, temp_all_covering_points_in_the_block;
	temp_all_S_in_the_block.clear(); temp_all_covering_points_in_the_block.clear();
	temp_all_S_in_the_block.resize(layer_graph.data.total_node_num);
	temp_all_covering_points_in_the_block.resize(layer_graph.data.total_node_num);

	std::cout << "ori_num_points_of_ori_in_all_the_area_S.size()" << ori_num_points_of_ori_in_all_the_area_S.size() << std::endl;

	std::vector<bool> area_s_layer_candidate(ori_num_points_of_ori_in_all_the_area_S.size(), false);
	for (int s_id = 0; s_id < static_cast<int>(ori_num_points_of_ori_in_all_the_area_S.size()); ++s_id) {
		const auto map_it = map_S_and_vertex.find(s_id);
		if (map_it == map_S_and_vertex.end() || judge_S_be_searched[s_id]) {
			continue;
		}

		bool is_currently_accessible = false;
		for (int count : ori_num_points_of_ori_in_all_the_area_S[s_id]) {
			if (count < 0) {
				cout << "error!!";
			}
			if (count == 0) {
				is_currently_accessible = true;
				break;
			}
		}
		if (is_currently_accessible) {
			cout << "Area S " << s_id << " is currently accessible." << std::endl;
			continue;
		}
		area_s_layer_candidate[s_id] = true;

		const int vertex_id = map_it->second;
		if (vertex_id < 0 || vertex_id >= static_cast<int>(V_2.size())) {
			continue;
		}
		const Eigen::Vector3d point(
			V_2[vertex_id](0, 0),
			V_2[vertex_id](1, 0),
			V_2[vertex_id](2, 0));
		const std::vector<int> layer_ids = FindContainingExtrudedMaterialLayers(
			layer_graph,
			vec_polygon,
			vec_hole_polygons,
			point,
			dh,
			layer_z_slack,
			layer_polygon_offset);
		for (int layer_id : layer_ids) {
			temp_all_S_in_the_block[layer_id].push_back(s_id);
		}
		if (!layer_ids.empty()) {
			judge_S_be_searched[s_id] = true;
		}
	}

	for (int covering_point_id = 0;
		covering_point_id < static_cast<int>(ori_all_the_covering_points.size());
		++covering_point_id) {
		if (covering_point_id >= static_cast<int>(judge_covering_points_be_searched.size())
			|| judge_covering_points_be_searched[covering_point_id]) {
			continue;
		}
		const auto map_it = map_covering_points_and_vertex.find(covering_point_id);
		if (map_it == map_covering_points_and_vertex.end()) {
			continue;
		}
		const int vertex_id = map_it->second;
		if (vertex_id < 0 || vertex_id >= static_cast<int>(V_2.size())) {
			continue;
		}
		const Eigen::Vector3d point(
			V_2[vertex_id](0, 0),
			V_2[vertex_id](1, 0),
			V_2[vertex_id](2, 0));
		const int layer_id = FindContainingExtrudedMaterialLayer(
			layer_graph,
			vec_polygon,
			vec_hole_polygons,
			point,
			dh,
			layer_z_slack,
			layer_polygon_offset);
		if (layer_id >= 0) {
			temp_all_covering_points_in_the_block[layer_id].push_back(covering_point_id);
			judge_covering_points_be_searched[covering_point_id] = true;
		}
	}



	int temp_all_S_in_the_block_count = 0;
	std::vector<Eigen::Vector3d> temp_all_S_in_the_block_points;
	std::vector<Eigen::Vector3d> assigned_area_s_layer_debug_points;
	std::vector<Eigen::Vector3d> unassigned_area_s_layer_debug_points;
	for (int j = 0; j < temp_all_S_in_the_block.size(); j++) {
		temp_all_S_in_the_block_count += temp_all_S_in_the_block[j].size();
		for (int s_id : temp_all_S_in_the_block[j]) {
			const auto map_it = map_S_and_vertex.find(s_id);
			if (map_it == map_S_and_vertex.end()) {
				continue;
			}
			const int vertex_id = map_it->second;
			if (vertex_id < 0 || vertex_id >= V.rows()) {
				continue;
			}
			temp_all_S_in_the_block_points.emplace_back(
				V(vertex_id, 0),
				V(vertex_id, 1),
				V(vertex_id, 2));
			if (vertex_id < V_2.size()) {
				assigned_area_s_layer_debug_points.emplace_back(
					V_2[vertex_id](0, 0),
					V_2[vertex_id](1, 0),
					V_2[vertex_id](2, 0));
			}
		}
	}
	for (int s_id = 0; s_id < static_cast<int>(area_s_layer_candidate.size()); ++s_id) {
		if (!area_s_layer_candidate[s_id] || judge_S_be_searched[s_id]) {
			continue;
		}
		const auto map_it = map_S_and_vertex.find(s_id);
		if (map_it == map_S_and_vertex.end()) {
			continue;
		}
		const int vertex_id = map_it->second;
		if (vertex_id < 0 || vertex_id >= static_cast<int>(V_2.size())) {
			continue;
		}
		unassigned_area_s_layer_debug_points.emplace_back(
			V_2[vertex_id](0, 0),
			V_2[vertex_id](1, 0),
			V_2[vertex_id](2, 0));
	}
	std::cout << "temp_all_S_in_the_block_count: " << temp_all_S_in_the_block_count << std::endl;
	//WriteDfsLayerContainmentDebugObj(
	//	UnaccessibleVisPath("access_debug_dfs_layer_area_s_relation_" + MakeAccessibilityDebugFileToken(file_name) + ".obj"),
	//	vec_polygon,
	//	vec_polygon_z,
	//	assigned_area_s_layer_debug_points,
	//	unassigned_area_s_layer_debug_points,
	//	std::max(0.12, dh * 0.08),
	//	layer_polygon_offset,
	//	layer_z_slack);
	//WriteDebugMarkersObj(
	//	UnaccessibleVisPath("access_debug_temp_all_S_in_the_block_" + MakeAccessibilityDebugFileToken(file_name) + ".obj"),
	//	temp_all_S_in_the_block_points,
	//	{ 0.05, 1.0, 1.0 },
	//	std::max(0.12, dh * 0.08));

	////////////////////////////////////////
	//std::cout << "DD" << endl;
	while (last_step_nodes.size() != 0 || candidate_nodes.size() != 0) {
		flag_continue = false;
		int now_last_node = last_step_nodes.front();
		bool jud_terminate = true;

		//update//
		while (tree_entries[now_last_node].parent != -1) {
			layer_graph.UpdateDegree_2(tree_entries[now_last_node].layer, -1);
			layer_graph.node_visited[tree_entries[now_last_node].layer] = true;
			now_last_node = tree_entries[now_last_node].parent;
		}
		//////////
		now_last_node = last_step_nodes.front();

		for (int i = 0; i < layer_graph.total_node_num; i++) {
			int v = i;
			bool flag_self_support = true;
			if (tree_entries.size() != 0 && tree_entries[now_last_node].layer == v) continue;
			bool jud_continue_2 = false;
			/*for (int j = 0; j < candidate_nodes.size(); j++)
				if (Tree_nodes[candidate_nodes[j]] == v) {
					jud_continue_2 = true;
					break;
				}*/
			if (jud_continue_2 == true) continue;
			if (layer_graph.node_visited[v]) continue;
			if (layer_graph.out_degree[v] != 0) continue;   //out-degree == 0
			//////////collision dependency edges///////////////
			if (tree_entries.size() != 0 && layer_graph.IsDepend_collision(tree_entries[now_last_node].layer, v))
				continue;
			///////////////////////////////////////////////////
			if (layer_graph.is_the_layer_self_suppot[v] == false && now_last_node != 0) {  //self-support constraint
				flag_self_support = false;
			}

			vector<int> all_S_in_the_block, all_covering_points_in_the_block;
			all_S_in_the_block = temp_all_S_in_the_block[i];
			all_covering_points_in_the_block = temp_all_covering_points_in_the_block[i];

			//try merge//
			bool jud_merge = true;
			for (int j = 0; j < all_S_in_the_block.size(); j++) {
				//check whether all S in the block is accessible//
				bool jud_merge_2 = false;
				for (int k = 0; k < ori_num_points_of_ori_in_all_the_area_S[all_S_in_the_block[j]].size(); k++) {
					if (ori_num_points_of_ori_in_all_the_area_S[all_S_in_the_block[j]][k] < 0)
						cout << "error!!";
					if (ori_num_points_of_ori_in_all_the_area_S[all_S_in_the_block[j]][k] == 0) {
						jud_merge_2 = true;
						break;
					}
				}
				if (jud_merge_2 == false)
					jud_merge = false;
				//////////////////////////////////////////////////
			}
			if (jud_merge == true && flag_self_support == false && previous_is_continue == false) {
				if (open_change_orientation == true)
					flag_continue = true;
			}
			if (jud_merge == true && flag_self_support == true) {
				jud_terminate = false;
				index_node++;
				DfsTreeEntry entry;
				entry.layer = v;
				entry.area_s = all_S_in_the_block;
				entry.sample_points = sample_point_in_layer[v];
				entry.covering_points = all_covering_points_in_the_block;
				entry.parent = now_last_node;
				tree_entries.push_back(entry);
				candidate_nodes.push_back(index_node);

				break; //DFS
			}
			////////////
		}
		if (jud_terminate == true) {
			if (now_last_node >= 0
				&& now_last_node < static_cast<int>(tree_entries.size())
				&& tree_entries[now_last_node].layer >= 0
				&& tree_entries[now_last_node].layer < layer_graph.total_node_num) {
				terminate_nodes.push_back(now_last_node);
			}
			else {
				std::cout << "[HybridManufacturing::DFS_search] No valid layer path remains; "
					<< "skip the root-only path." << std::endl;
			}
			break;
		}

		//update//
		while (tree_entries[now_last_node].parent != -1) {
			layer_graph.UpdateDegree_2(tree_entries[now_last_node].layer, 1);
			layer_graph.node_visited[tree_entries[now_last_node].layer] = false;
			now_last_node = tree_entries[now_last_node].parent;
		}
		//////////

		//sort_candidate_nodes(candidate_nodes, Tree_nodes_for_S);
		last_step_nodes.pop();
		//cout << candidate_nodes.size() << endl;
		if (last_step_nodes.size() == 0) {
			int cont_w = 0;
			while (candidate_nodes.size() != 0 && cont_w < W2 && cont_w < candidate_nodes.size()) {
				last_step_nodes.push(candidate_nodes[cont_w]);
				cont_w++;
			}
			candidate_nodes.clear();
		}
	}

	//every path (blcok) as a candidate node of outer-beam search (save slice layers)//
	//final_pathes save the index of slice layer (layer_graph.total_node_num)
	pathes_include_S.clear();
	pathes_include_sample_points.clear();
	paths_include_covering_points.clear();
	vector<vector<int>> final_pathes;
	vector<int> final_pathes_height;
	for (int i = 0; i < terminate_nodes.size(); i++) {
		if (terminate_nodes[i] < 0
			|| terminate_nodes[i] >= static_cast<int>(tree_entries.size())
			|| tree_entries[terminate_nodes[i]].layer < 0
			|| tree_entries[terminate_nodes[i]].layer >= layer_graph.total_node_num) {
			continue;
		}
		int height = 0;
		vector<int> current_path, current_path_include_S, current_path_include_sample_points, current_path_include_covering_points;
		current_path.clear();
		current_path_include_S.clear();
		current_path_include_sample_points.clear();
		current_path_include_covering_points.clear();
		int current_node = terminate_nodes[i];
		current_path.push_back(tree_entries[terminate_nodes[i]].layer);
		for (int j = 0; j < tree_entries[terminate_nodes[i]].area_s.size(); j++)
			current_path_include_S.push_back(tree_entries[terminate_nodes[i]].area_s[j]);
		for (int j = 0; j < tree_entries[terminate_nodes[i]].sample_points.size(); j++)
			current_path_include_sample_points.push_back(tree_entries[terminate_nodes[i]].sample_points[j]);
		for (int j = 0; j < tree_entries[terminate_nodes[i]].covering_points.size(); j++)
			current_path_include_covering_points.push_back(tree_entries[terminate_nodes[i]].covering_points[j]);
		while (tree_entries[current_node].parent > 0) {
			current_node = tree_entries[current_node].parent;
			current_path.push_back(tree_entries[current_node].layer);
			height++;
			for (int j = 0; j < tree_entries[current_node].area_s.size(); j++) {
				current_path_include_S.push_back(tree_entries[current_node].area_s[j]);
			}

			for (int j = 0; j < tree_entries[current_node].sample_points.size(); j++)
				current_path_include_sample_points.push_back(tree_entries[current_node].sample_points[j]);
			for (int j = 0; j < tree_entries[current_node].covering_points.size(); j++)
				current_path_include_covering_points.push_back(tree_entries[current_node].covering_points[j]);
		}
		std::reverse(current_path.begin(), current_path.end());
		std::reverse(current_path_include_S.begin(), current_path_include_S.end());
		std::reverse(current_path_include_sample_points.begin(), current_path_include_sample_points.end());
		std::reverse(current_path_include_covering_points.begin(), current_path_include_covering_points.end());
		final_pathes.push_back(current_path);
		final_pathes_height.push_back(height);
		pathes_include_S.push_back(current_path_include_S);
		pathes_include_sample_points.push_back(current_path_include_sample_points);
		paths_include_covering_points.push_back(current_path_include_covering_points);
	}

	//sort final_pathes and select W2 pathes//
	for (int i = 0; i < final_pathes.size(); i++)
		for (int j = i + 1; j < final_pathes.size(); j++) {
			if (final_pathes_height[i] < final_pathes_height[j]) {
				swap(final_pathes[i], final_pathes[j]);
				swap(final_pathes_height[i], final_pathes_height[j]);
				swap(pathes_include_S[i], pathes_include_S[j]);
				swap(pathes_include_sample_points[i], pathes_include_sample_points[j]);
				swap(paths_include_covering_points[i], paths_include_covering_points[j]);
			}
		}
	while (final_pathes.size() > W2) {
		auto itr = final_pathes.begin() + W2;
		final_pathes.erase(itr);
		auto itr2 = pathes_include_S.begin() + W2;
		pathes_include_S.erase(itr2);
		auto itr3 = pathes_include_sample_points.begin() + W2;
		pathes_include_sample_points.erase(itr3);
		auto itr4 = paths_include_covering_points.begin() + W2;
		paths_include_covering_points.erase(itr4);
	}
	////////////////////

	all_cut_layers = FindAllCutLayers(layer_graph, final_pathes, all_cut_layers_dependency_layer, jud_admit);

	vector<vector<vector<cv::Point3d>>> all_selected_points;	//存储所有被选择的路径上的切片轮廓的各个点的坐标
	vector<vector<vector<cv::Point3d>>> all_selected_points_contain;
	vector<CutLayerHoles> all_selected_layer_holes;
	for (int i = 0; i < final_pathes.size(); i++) {
		vector<vector<cv::Point3d>> temp_vec_1;
		all_selected_points.push_back(temp_vec_1);
		all_selected_points_contain.push_back(temp_vec_1);
		all_selected_layer_holes.emplace_back();
		for (int j = 0; j < final_pathes[i].size(); j++) {
			vector<cv::Point3d> temp_vec_2;
			all_selected_points[i].push_back(temp_vec_2);
			all_selected_points_contain[i].push_back(temp_vec_2);
			all_selected_layer_holes[i].emplace_back();
			const int layer_id = final_pathes[i][j];
			const auto index_it = layer_graph.data.index.find(layer_id);
			if (layer_id < 0
				|| layer_id >= layer_graph.total_node_num
				|| index_it == layer_graph.data.index.end()) {
				std::cout << "[HybridManufacturing::DFS_search] Invalid layer id in final path: "
					<< layer_id << std::endl;
				jud_admit = false;
				return;
			}
			const pair<int, int> index_slice_point = index_it->second;
			if (index_slice_point.first < 0
				|| index_slice_point.first >= static_cast<int>(layer_graph.data.slice_points.size())
				|| index_slice_point.first >= static_cast<int>(layer_graph.data.slice_points_contain.size())
				|| index_slice_point.first >= static_cast<int>(layer_graph.data.slice_points_holes.size())
				|| index_slice_point.first >= static_cast<int>(layer_graph.data.z_value.size())
				|| index_slice_point.second < 0
				|| index_slice_point.second >= static_cast<int>(layer_graph.data.slice_points[index_slice_point.first].size())
				|| index_slice_point.second >= static_cast<int>(layer_graph.data.slice_points_contain[index_slice_point.first].size())
				|| index_slice_point.second >= static_cast<int>(layer_graph.data.slice_points_holes[index_slice_point.first].size())
				|| index_slice_point.second >= static_cast<int>(layer_graph.data.z_value[index_slice_point.first].size())
				|| layer_graph.data.z_value[index_slice_point.first][index_slice_point.second].empty()) {
				std::cout << "[HybridManufacturing::DFS_search] Invalid slice index for layer "
					<< layer_id
					<< ": (" << index_slice_point.first
					<< ", " << index_slice_point.second << ")"
					<< std::endl;
				jud_admit = false;
				return;
			}
			const auto& selected_slice_points =
				layer_graph.data.slice_points[index_slice_point.first][index_slice_point.second];
			const auto& selected_z_values =
				layer_graph.data.z_value[index_slice_point.first][index_slice_point.second];
			const std::size_t selected_point_count = std::min(
				selected_slice_points.size(),
				selected_z_values.size());
			if (selected_point_count != selected_slice_points.size()) {
				std::cout << "[HybridManufacturing::DFS_search] Slice point/z-value count mismatch for layer "
					<< layer_id
					<< ": points=" << selected_slice_points.size()
					<< ", z_values=" << selected_z_values.size()
					<< std::endl;
			}
			for (std::size_t k = 0; k < selected_point_count; ++k) {
				cv::Point3d current_point(
					selected_slice_points[k].x(),
					selected_slice_points[k].y(),
					selected_z_values[k]);
				all_selected_points[i][j].push_back(current_point);
			}
			const double selected_z = selected_z_values.front();
			const auto& selected_holes =
				layer_graph.data.slice_points_holes[index_slice_point.first][index_slice_point.second];
			for (const auto& hole : selected_holes) {
				all_selected_layer_holes[i][j].emplace_back();
				auto& output_hole = all_selected_layer_holes[i][j].back();
				output_hole.reserve(hole.size());
				for (const Eigen::Vector2d& point : hole) {
					output_hole.emplace_back(point.x(), point.y(), selected_z);
				}
			}
			if (!all_selected_layer_holes[i][j].empty()) {
				all_selected_points_contain[i][j] = all_selected_layer_holes[i][j].front();
			}
		}
	}


	//////////visualization all triangles cross by layers//////////
	//vector<vector<vector<Vertex*>>> show_triangles(final_pathes[0].size());
	//for (int i = 0; i < 1; i++) {
	//	for (int j = 0; j < final_pathes[0].size(); j++) {
	//		show_triangles[j].resize(layer_graph.all_triangles_of_layers[final_pathes[i][j]].size());
	//		for (int k = 0; k < layer_graph.all_triangles_of_layers[final_pathes[i][j]].size(); k++) {
	//			Vertex* v1 = layer_graph.all_triangles_of_layers[final_pathes[i][j]][k]->vertices_2[0];
	//			Vertex* v2 = layer_graph.all_triangles_of_layers[final_pathes[i][j]][k]->vertices_2[1];
	//			Vertex* v3 = layer_graph.all_triangles_of_layers[final_pathes[i][j]][k]->vertices_2[2];
	//			//double ans = (v2->x - v1->x) * (v2->y - v3->y) - (v2->y - v1->y) * (v2->x - v3->x);
	//			//if (ans > 0)	//is clockwise
	//			//	swap(v2, v3);
	//			show_triangles[j][k].push_back(v1);
	//			show_triangles[j][k].push_back(v2);
	//			show_triangles[j][k].push_back(v3);
	//		}
	//	}
	//}
	//std::ofstream dstream(".\\vis\\show_triangles.stl");
	//if (!dstream.is_open()) {
	//	std::cout << "can not open " << std::endl;
	//}
	//dstream << "solid STL generated by MeshLab" << std::endl;
	//for (int j = 0; j < final_pathes[0].size(); j++) {
	//	for (int k = 0; k < show_triangles[j].size(); k++) {
	//		dstream << "  facet normal " << "0 0 0" << std::endl;
	//		dstream << "    outer loop" << std::endl;
	//		for (int l = 0; l < 3; l++) {
	//			dstream << "      vertex  " << show_triangles[j][k][l]->x << " " << show_triangles[j][k][l]->y << " " << show_triangles[j][k][l]->z << std::endl;
	//		}
	//		dstream << "    endloop" << std::endl;
	//		dstream << "  endfacet" << std::endl;
	//	}
	//}
	//dstream << "endsolid vcg" << std::endl;
	//dstream.close();

	all_solutions_of_selected_layers = all_selected_points;
	all_solutions_of_selected_layers_contain = all_selected_points_contain;
	all_solutions_of_selected_layer_holes = all_selected_layer_holes;
	return;
	//////////////////////////////////////////////////////////////////////////////////
}

void HybridManufacturing::sort_candidate_nodes(
	vector<int>& candidate_nodes,
	const vector<BeamTreeEntry>& tree_entries,
	vector<vector<int>> final_pathes_include_S,
	OrientationScores& all_calculated_value,
	vector<vector<int>> final_pathes_include_covering_points,
	int height_of_beam_search,
	vector<vector<Eigen::Vector3d>> save_ori,
	OrientationScores& pure_value,
	int id_continue)
{
	double W_less_area_S = 0, W_more_slice_layers = 0.4, W_covering_points = 0.6, W_less_clipping_plane = 0, W_fragile = 0, W_larger_base = 0, W_orientation = 0, W_projected = 0;

	//double W_less_area_S = 0, W_more_slice_layers = 0.1, W_covering_points = 0.1, W_less_clipping_plane = 0.4, W_fragile = 0, W_larger_base = 0, W_orientation = 0.4, W_projected = 0;
	//double W_less_area_S = 0, W_more_slice_layers = 0, W_covering_points = 0, W_less_clipping_plane = 0, W_fragile = 0, W_larger_base = 0, W_orientation = 0, W_projected = 1;

	//if(height_of_beam_search >= 12)  //test for figure of paper
	//	W_less_area_S = 0.1, W_more_slice_layers = 1, W_covering_points = 0.1, W_less_clipping_plane = 0, W_fragile = 0, W_larger_base = 0, W_orientation = 0, W_projected = 0;
	/*if(height_of_beam_search == 2)
		W_less_area_S = 0.9, W_more_slice_layers = 0.1, W_covering_points = 0, W_less_clipping_plane = 0, W_fragile = 0, W_larger_base = 0, W_orientation = 0, W_projected = 0;*/

		//random test
		//int rand_id = rand() % all_calculated_value.size();
		/*int rand_id = 5;
		swap(candidate_nodes[0], candidate_nodes[rand_id]);
		swap(all_calculated_value[0], all_calculated_value[rand_id]);
		return;*/

	// for (int i = 0; i < all_calculated_value.size(); i++) {
	// 	if (all_calculated_value[i].number_of_remaining_face < terminate_threshold_of_number_of_faces && save_ori[candidate_nodes[i]][0].x() == 0 && save_ori[candidate_nodes[i]][0].y() == 0) {
	// 		swap(candidate_nodes[0], candidate_nodes[i]);
	// 		swap(all_calculated_value[0], all_calculated_value[i]);
	// 		cout << "Remaining face: " << all_calculated_value[i].number_of_remaining_face << endl;
	// 		return;
	// 	}
	// }
	// for (int i = 0; i < all_calculated_value.size(); i++) {
	// 	if (height_of_beam_search == 4 && save_ori[candidate_nodes[i]][0].x() == 0 && save_ori[candidate_nodes[i]][0].y() == 0) {
	// 		swap(candidate_nodes[0], candidate_nodes[i]);
	// 		swap(all_calculated_value[0], all_calculated_value[i]);
	// 		cout << "Remaining face: " << all_calculated_value[i].number_of_remaining_face << endl;
	// 		return;
	// 	}
	// }

	/*if(height_of_beam_search == 2)
		W_less_area_S = 0, W_more_slice_layers = 0.1, W_covering_points = 0.1, W_less_clipping_plane = 0.8, W_fragile = 0, W_larger_base = 0;*/

		//normalization - local//
		/*int max_area_S = -MAX_I, min_area_S = MAX_I;
		int max_slice_layers = -MAX_I, min_slice_layers = MAX_I;
		int max_clipping_plane = -MAX_I, min_clipping_plane = MAX_I;
		double max_large_base = -MAX_I, min_large_base = MAX_I;
		for (int i = 0; i < all_calculated_value.size(); i++) {
			max_large_base = max(max_large_base, all_calculated_value[i].large_base);
			min_large_base = min(min_large_base, all_calculated_value[i].large_base);

			max_area_S = max(max_area_S, int(final_pathes_include_S[candidate_nodes[i]].size()));
			min_area_S = min(min_area_S, int(final_pathes_include_S[candidate_nodes[i]].size()));

			max_slice_layers = max(max_slice_layers, int(tree_entries[candidate_nodes[i]].layers.size()));
			min_slice_layers = min(min_slice_layers, int(tree_entries[candidate_nodes[i]].layers.size()));

			max_clipping_plane = max(max_clipping_plane, int(tree_entries[candidate_nodes[i]].cut_layers.size()));
		}


		for (int i = 0; i < all_calculated_value.size(); i++) {
			if (max_large_base - min_large_base != 0)
				all_calculated_value[i].large_base = (all_calculated_value[i].large_base - min_large_base) / (max_large_base - min_large_base);
			else
				all_calculated_value[i].large_base = 0;

			if (max_area_S - min_area_S != 0)
				all_calculated_value[i].value_of_area_S = double(final_pathes_include_S[candidate_nodes[i]].size() - min_area_S) / double(max_area_S - min_area_S);
			else
				all_calculated_value[i].value_of_area_S = 0;

			if (max_slice_layers - min_slice_layers != 0)
			all_calculated_value[i].value_of_more_slice_layers = double(tree_entries[candidate_nodes[i]].layers.size() - min_slice_layers) / double(max_slice_layers - min_slice_layers);
			else
				all_calculated_value[i].value_of_more_slice_layers = 0;

			if (max_clipping_plane - min_clipping_plane != 0)
			all_calculated_value[i].value_of_less_clipping_plane = 1 - double(tree_entries[candidate_nodes[i]].cut_layers.size() - min_clipping_plane) / double(max_clipping_plane - min_clipping_plane);
			else
				all_calculated_value[i].value_of_less_clipping_plane = 0;
		}*/
		/////////////////


		//Calculate the accumulated value//
	int* sum_area_S = new int[all_calculated_value.size()];
	double* sum_slice_layers = new double[all_calculated_value.size()];
	int* sum_clipping_plane = new int[all_calculated_value.size()];
	double* sum_larger_base = new double[all_calculated_value.size()];
	double* sum_covering_points = new double[all_calculated_value.size()];
	double* sum_fragile = new double[all_calculated_value.size()];
	double* sum_orientation = new double[all_calculated_value.size()];
	double* sum_projected = new double[all_calculated_value.size()];
	for (int i = 0; i < all_calculated_value.size(); i++) {
		sum_area_S[i] = 0;
		sum_slice_layers[i] = 0;
		sum_clipping_plane[i] = 0;
		sum_larger_base[i] = 0;
		sum_covering_points[i] = 0;
		sum_orientation[i] = 0;
		sum_fragile[i] = all_calculated_value[i].value_of_fragile;
		//cout << "((((((" << sum_fragile[i] << endl;
		sum_projected[i] = all_calculated_value[i].value_of_projected;
		int index_node = candidate_nodes[i];
		const int parent_id = tree_entries[index_node].parent_id;
		if (parent_id >= 0) {
			pure_value[i].value_of_orientation = save_ori[index_node][0].dot(save_ori[parent_id][0]);
		}
		else {
			pure_value[i].value_of_orientation = 0;
		}
		pure_value[i].value_of_fragile = all_calculated_value[i].value_of_fragile;
		pure_value[i].value_of_projected = all_calculated_value[i].value_of_projected;
		while (tree_entries[index_node].parent_id != -1) {
			sum_area_S[i] += int(final_pathes_include_S[index_node].size());
			for (int j = 0; j < tree_entries[index_node].layers.size(); j++) {
				for (int k = 0; k < tree_entries[index_node].layers[j].size() - 1; k++) {
					sum_slice_layers[i] += distance3d(tree_entries[index_node].layers[j][k], tree_entries[index_node].layers[j][k + 1]);
				}
			}
			sum_clipping_plane[i] += int(tree_entries[index_node].cut_layers.size());
			sum_larger_base[i] += tree_entries[index_node].larger_base;

			//sum_covering_points[i] += int(final_pathes_include_covering_points[index_node].size());
			for (int j = 0; j < final_pathes_include_covering_points[index_node].size(); j++)
				sum_covering_points[i] += int(all_the_covering_points[final_pathes_include_covering_points[index_node][j]].size());
			index_node = tree_entries[index_node].parent_id;
		}
		index_node = candidate_nodes[i];
		while (tree_entries[index_node].parent_id > 0) {
			sum_orientation[i] += save_ori[index_node][0].dot(save_ori[tree_entries[index_node].parent_id][0]);
			index_node = tree_entries[index_node].parent_id;
		}
	}
	////////////////////////////////////


	//normalization - accumulated//
	int max_area_S = -MAX_I, min_area_S = MAX_I;
	double max_slice_layers = -MAX_I, min_slice_layers = MAX_I;
	int max_clipping_plane = -MAX_I, min_clipping_plane = MAX_I;
	double max_large_base = -MAX_I, min_large_base = MAX_I;
	double max_covering_points = -MAX_I, min_covering_points = MAX_I;
	double max_fragile = -MAX_I, min_fragile = MAX_I;
	double max_orientation = -MAX_I, min_orientation = MAX_I;
	double max_projected = -MAX_I, min_projected = MAX_I;
	for (int i = 0; i < all_calculated_value.size(); i++) {
		max_large_base = max(max_large_base, sum_larger_base[i]);
		min_large_base = min(min_large_base, sum_larger_base[i]);

		max_area_S = max(max_area_S, sum_area_S[i]);
		min_area_S = min(min_area_S, sum_area_S[i]);

		max_slice_layers = max(max_slice_layers, sum_slice_layers[i]);
		min_slice_layers = min(min_slice_layers, sum_slice_layers[i]);

		max_clipping_plane = max(max_clipping_plane, sum_clipping_plane[i]);
		min_clipping_plane = min(min_clipping_plane, sum_clipping_plane[i]);

		max_covering_points = max(max_covering_points, sum_covering_points[i]);
		min_covering_points = min(min_covering_points, sum_covering_points[i]);

		max_fragile = max(max_fragile, sum_fragile[i]);
		min_fragile = min(min_fragile, sum_fragile[i]);

		max_orientation = max(max_orientation, sum_orientation[i]);
		min_orientation = min(min_orientation, sum_orientation[i]);

		max_projected = max(max_projected, sum_projected[i]);
		min_projected = min(min_projected, sum_projected[i]);
	}

	for (int i = 0; i < all_calculated_value.size(); i++) {
		if (max_large_base - min_large_base != 0)
			all_calculated_value[i].large_base = (sum_larger_base[i] - min_large_base) / (max_large_base - min_large_base);
		else
			all_calculated_value[i].large_base = 1;  //全部暂时改为1，之后还是需要-1用于标识

		if (max_area_S - min_area_S != 0)
			all_calculated_value[i].value_of_area_S = 1 - double(sum_area_S[i] - min_area_S) / double(max_area_S - min_area_S);
		else
			all_calculated_value[i].value_of_area_S = 1;

		if (max_slice_layers - min_slice_layers != 0)
			all_calculated_value[i].value_of_more_slice_layers = double(sum_slice_layers[i] - min_slice_layers) / double(max_slice_layers - min_slice_layers);
		else
			all_calculated_value[i].value_of_more_slice_layers = 1;

		if (max_clipping_plane - min_clipping_plane != 0)
			all_calculated_value[i].value_of_less_clipping_plane = 1 - double(sum_clipping_plane[i] - min_clipping_plane) / double(max_clipping_plane - min_clipping_plane);
		else
			all_calculated_value[i].value_of_less_clipping_plane = 1;

		if (max_covering_points - min_covering_points != 0)
			all_calculated_value[i].value_of_covering_points = double(sum_covering_points[i] - min_covering_points) / double(max_covering_points - min_covering_points);
		else
			all_calculated_value[i].value_of_covering_points = 1;

		if (max_fragile - min_fragile != 0)
			all_calculated_value[i].value_of_fragile = 1 - double(sum_fragile[i] - min_fragile) / double(max_fragile - min_fragile);
		else
			all_calculated_value[i].value_of_fragile = 1;

		if (max_orientation - min_orientation != 0)
			all_calculated_value[i].value_of_orientation = double(sum_orientation[i] - min_orientation) / double(max_orientation - min_orientation);
		else
			all_calculated_value[i].value_of_orientation = 1;

		/*if (height_of_beam_search == 3) {
			if (max_projected - min_projected != 0)
				all_calculated_value[i].value_of_projected = double(sum_projected[i] - min_projected) / double(max_projected - min_projected);
			else
				all_calculated_value[i].value_of_projected = -1;
		}*/
		//else {
		if (max_projected - min_projected != 0)
			all_calculated_value[i].value_of_projected = 1 - double(sum_projected[i] - min_projected) / double(max_projected - min_projected);
		else
			all_calculated_value[i].value_of_projected = 1;
		//}
	}
	/////////////////

	//sort
	for (int i = 0; i < candidate_nodes.size(); i++) {
		for (int j = i + 1; j < candidate_nodes.size(); j++) {
			auto calc_score = [&](int index) {
				const auto& value = all_calculated_value[index];
				return value.value_of_area_S * W_less_area_S
					+ value.value_of_more_slice_layers * W_more_slice_layers
					+ value.value_of_less_clipping_plane * W_less_clipping_plane
					+ value.large_base * W_larger_base
					+ value.value_of_covering_points * W_covering_points
					+ value.value_of_fragile * W_fragile
					+ value.value_of_orientation * W_orientation
					+ value.value_of_projected * W_projected
					+ value.value_of_sub_patches * 0.0;
				};

			if (calc_score(i) < calc_score(j)) {
				swap(candidate_nodes[i], candidate_nodes[j]);
				swap(all_calculated_value[i], all_calculated_value[j]);
				swap(pure_value[i], pure_value[j]);
			}
		}
	}


	//if (height_of_beam_search == 6) {
		//cout << endl << "&&" << all_calculated_value[0].number_of_remaining_face << endl;
		//cout << endl << "&&" << all_calculated_value[1].number_of_remaining_face << endl;
		//}

	//先删除
	cout << "#####################" << height_of_beam_search << " " << id_continue << "#####################" << endl;
	//if (height_of_beam_search == 2) {  //julia_vase
	//	int cont = 0;
	//	while (cont < 2) {  //5
	//		candidate_nodes.erase(candidate_nodes.begin());
	//		all_calculated_value.erase(all_calculated_value.begin());
	//		cont++;
	//	}
	//}
	// 
	//if (height_of_beam_search == 3) {
	//	int cont = 0;
	//	while (cont < 1) {  //5
	//		candidate_nodes.erase(candidate_nodes.begin());
	//		all_calculated_value.erase(all_calculated_value.begin());
	//		cont++;
	//	}
	//}
	//else if (height_of_beam_search == 5) {
	//	int cont = 0;
	//	while (cont <1) {  //6
	//		candidate_nodes.erase(candidate_nodes.begin());
	//		all_calculated_value.erase(all_calculated_value.begin());
	//		cont++;
	//	}
	//}
	//if (height_of_beam_search == 2 && id_continue == 0) {
	//	int cont = 0;
	//	while (cont < 1) {  //6
	//		candidate_nodes.erase(candidate_nodes.begin());
	//		all_calculated_value.erase(all_calculated_value.begin());
	//		cont++;
	//	}
	//}

	//if (height_of_beam_search == 4) {
	//	int cont = 0;
	//	while (cont < 5) {  //5
	//		candidate_nodes.erase(candidate_nodes.begin());
	//		all_calculated_value.erase(all_calculated_value.begin());
	//		cont++;
	//	}
	//}
	//if (height_of_beam_search == 5) {
	//	int cont = 0;
	//	while (cont < 8) {  //5
	//		candidate_nodes.erase(candidate_nodes.begin());
	//		all_calculated_value.erase(all_calculated_value.begin());
	//		cont++;
	//	}
	//}
	//if (height_of_beam_search == 6) {
	//	int cont = 0;
	//	while (cont < 13) {  //5
	//		candidate_nodes.erase(candidate_nodes.begin());
	//		all_calculated_value.erase(all_calculated_value.begin());
	//		cont++;
	//	}
	//}
	//if (height_of_beam_search == 7) {
	//	int cont = 0;
	//	while (cont < 0) {  //5
	//		candidate_nodes.erase(candidate_nodes.begin());
	//		all_calculated_value.erase(all_calculated_value.begin());
	//		cont++;
	//	}
	//}
	//if (height_of_beam_search == 13) {
	//	int cont = 0;
	//	while (cont < 1) {  //5
	//		candidate_nodes.erase(candidate_nodes.begin());
	//		all_calculated_value.erase(all_calculated_value.begin());
	//		cont++;
	//	}
	//}
	//if (height_of_beam_search == 14) {
	//	int cont = 0;
	//	while (cont < 2) {  //5
	//		candidate_nodes.erase(candidate_nodes.begin());
	//		all_calculated_value.erase(all_calculated_value.begin());
	//		cont++;
	//	}
	//}
	//if (height_of_beam_search == 18) {
	//	int cont = 0;
	//	while (cont < 1) {  //5
	//		candidate_nodes.erase(candidate_nodes.begin());
	//		all_calculated_value.erase(all_calculated_value.begin());
	//		cont++;
	//	}
	//}
}

void HybridManufacturing::sort_candidate_nodes(vector<int>& candidate_nodes, vector<vector<int>> Tree_nodes_for_S)
{
	/*for (int i = 0; i < candidate_nodes.size(); i++) {
		for (int j = i + 1; j < candidate_nodes.size(); j++) {
			if (Tree_nodes_for_S[candidate_nodes[i]].size() < Tree_nodes_for_S[candidate_nodes[j]].size()) {
				swap(candidate_nodes[i], candidate_nodes[j]);
			}
		}
	}*/
}

void HybridManufacturing::subtractive_remove_output(const vector<TRiangle>& need_detect_triangle, const Slicer_2& current_slicer, int height_of_beam_search, int cont_number_of_queue)
{
	std::string filename = PatchVisPath(
		"block_patch-" + to_string(height_of_beam_search)
		+ "_" + to_string(cont_number_of_queue)
		+ ".obj");
	std::ofstream file(filename);
	if (!file)
	{
		std::cout << "subtractive_remove_output !file" << std::endl;
	}

	for (auto& p : current_slicer.positions) {
		file << "v " << p[0] << " " << p[1] << " " << p[2] << "\n";
	}
	for (auto t : need_detect_triangle) {
		for (int i = 0; i < 3; ++i) t[i]++;
		file << "f " << t[0] << " " << t[1] << " " << t[2] << "\n";
	}
}

bool HybridManufacturing::CheckToolCollisionWithCell(
	const Eigen::Vector3d& center_point,
	const std::vector<Eigen::Vector3d>& target_cell_vertices,
	double max_z_target,
	const cutter& nozzle,
	double z_threshold_divisor,
	double xy_tolerance) const
{
	(void)xy_tolerance;

	const double z_threshold = nozzle.cylinder_r / z_threshold_divisor;
	const double height_diff = max_z_target - center_point.z();

	// 快速排除：目标在刀尖球下方
	if (height_diff <= z_threshold) {
		return false;
	}

	// 快速确认：目标超出刀具最大高度
	if (height_diff > nozzle.total_height) {
		return true;
	}

	// 使用第一个顶点进行粗筛
	if (!target_cell_vertices.empty()) {
		const double dx = target_cell_vertices[0](0, 0) - center_point.x();
		const double dy = target_cell_vertices[0](1, 0) - center_point.y();
		const double dist_xy_sq = dx * dx + dy * dy;

		if (height_diff > nozzle.cylinder_height_threshold) {
			if (dist_xy_sq > nozzle.carriage_check_radius_sq) {
				return false;
			}
		}
		else {
			if (dist_xy_sq > nozzle.cylinder_check_radius_sq) {
				return false;
			}
		}
	}

	// 精确检测每个边界顶点
	for (const auto& vertex : target_cell_vertices) {
		const double vx = vertex(0, 0);
		const double vy = vertex(1, 0);
		const double vz = vertex(2, 0);
		const double diff_z = vz - center_point.z();

		if (diff_z <= z_threshold) {
			continue;
		}

		const double dx = vx - center_point.x();
		const double dy = vy - center_point.y();
		const double dist_xy_sq = dx * dx + dy * dy;

		if (diff_z <= nozzle.cylinder_height_threshold) { //nozzle.cylinder_r + nozzle.cylinder_height
			if (dist_xy_sq < nozzle.cylinder_r_sq) {
				return true;
			}
		}
		else if (diff_z <= nozzle.total_height) { //nozzle.cylinder_r + nozzle.cylinder_height + nozzle.carriage_height
			if (dist_xy_sq < nozzle.carriage_r_sq) {
				return true;
			}
		}
	}

	return false;
}
