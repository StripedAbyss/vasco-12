#pragma once
#include <vector>
#include <array>
#include <string>
#include <map>
#include <regex>
#include <sstream>
#include <fstream>
#include <cmath>

namespace vasco {

class Slicer {
public:
    using Vec3    = std::array<double, 3>;
    using Tri3    = std::array<int, 3>;

    std::vector<Vec3>  positions;
    std::vector<Tri3>  triangles;
    std::vector<bool>  jud_plane;          // 标记由切割产生的新三角形

    Vec3 origin{0.0, 0.0, 0.0};
    Vec3 normal{0.0, 0.0, 1.0};

    void clear() {
        positions.clear();
        triangles.clear();
        jud_plane.clear();
    }

	bool load(const std::string& filename) {
		clear();
		std::ifstream file(filename);
		if (!file) return false;
		std::string line;
		while (std::getline(file, line)) {
			std::istringstream iss(line);
			std::string prefix;
			iss >> prefix;
			if (prefix == "v") {
				Vec3 p{};
				iss >> p[0] >> p[1] >> p[2];
				positions.push_back(p);
			}
			else if (prefix == "f") {
				Tri3 t{};
				iss >> t[0] >> t[1] >> t[2];
				for (int i = 0; i < 3; ++i)
					t[i]--;
				triangles.push_back(t);
			}
		}
		return true;
	}

    bool save(const std::string& filename) const {
        std::ofstream file(filename);
        if (!file) return false;
        for (auto& p : positions) {
            file << "v " << p[0] << " " << p[1] << " " << p[2] << "\n";
        }
        for (auto t : triangles) {
            for (int i = 0; i < 3; ++i) t[i]++;
            file << "f " << t[0] << " " << t[1] << " " << t[2] << "\n";
        }
        return true;
    }

    void cut() {
        // We loop on all triangle indices (tid) in the mesh.
        // The do_triangle function is called on each tid,
        // and will split the triangle when required.
        //
        // If do_triangle returns false, the triangle doesn't intersect the mesh.
        // The mesh (and triangle) is unchanged and we go the next triangle (tid+=1).
        //
        // If do_triangle return true, the triangle intersects the mesh.
        // The triangle has been modified and a new triangle has been appended to the mesh.
        // Since the modified triangle might still require splitting, we stay on it (tid+=0).
        // The new appended triangle will be processed at the end.
        //
        // See do_triangle.png for a picture.
        intersections.clear();
        if (origin[0] == 0 && origin[1] == 0 && origin[2] == 0) return;
        for (size_t tid = 0; tid < triangles.size(); tid += !splitTriangle(tid));
    }

    bool read_json(const std::string& filename) {
        //
        // Read JSON file with plane coordinates.
        //
        std::ifstream file(filename);
        if (!file) return false;
        std::stringstream buffer;
        buffer << file.rdbuf();
        std::string json = buffer.str();

        std::string o = std::regex_replace(json,
            std::regex(".*\"origin\"[^[]+\\[", std::regex::extended), "");
        o = std::regex_replace(o, std::regex("(\\].*)|,", std::regex::extended), "");
        std::istringstream(o) >> origin[0] >> origin[1] >> origin[2];

        std::string n = std::regex_replace(json,
            std::regex(".*\"normal\"[^[]+\\[", std::regex::extended), "");
        n = std::regex_replace(n, std::regex("(\\].*)|,", std::regex::extended), "");
        std::istringstream(n) >> normal[0] >> normal[1] >> normal[2];
        return true;
    }

private:
    // Floating point arithmetic is non-exact which might
    // cause problems when computing intersections.
    // We fix that using a hardcoded precision value.
    static constexpr double precision = 1e-7;

    // Keep track of indices of intersection points.
    std::map<std::pair<int,int>, int> intersections;

    //
    // Compute lambda such that lambda*P + (1-lambda)*Q
    // is the intersection between the plane and line PQ.
    //      lambda is in [0, 1]   =>  [PQ] intersects plane
    //      lambda is infinite    =>  [PQ] parallel to plane
    //      lambda is NaN         =>  [PQ] contained in plane
    //
    double computeLambda(const Vec3& p, const Vec3& q) const {
        double num = (origin[0] - q[0]) * normal[0]
                   + (origin[1] - q[1]) * normal[1]
                   + (origin[2] - q[2]) * normal[2];
        double den = (p[0] - q[0]) * normal[0]
                   + (p[1] - q[1]) * normal[1]
                   + (p[2] - q[2]) * normal[2];
        if (std::abs(num) < 1e-10 || std::abs(den) < 1e-10) return 0.0;
        return den == 0.0 ? 0.0 : num / den;
    }

    //
    // If edge [ij] intersects plane, add intersection vertex to mesh and return its index.
    // If the intersection vertex has already been computed before, only return its index.
    // If intersection does not exist, return -1.
    //
    int intersection(int i, int j) {
        if (i > j) std::swap(i, j);
        const auto key = std::make_pair(i, j);
        auto found = intersections.find(key);
        if (found != intersections.end()) return found->second;

        const Vec3& p = positions[i];
        const Vec3& q = positions[j];
        double lambda = computeLambda(p, q);
        if (!std::isfinite(lambda) || lambda < precision || lambda > 1 - precision) return -1;

        int m = static_cast<int>(positions.size());
        positions.push_back({
            lambda * p[0] + (1 - lambda) * q[0],
            lambda * p[1] + (1 - lambda) * q[1],
            lambda * p[2] + (1 - lambda) * q[2]
        });
        intersections[key] = m;
        return m;
    }

    //
    // If triangle tid intersects plane, split it and returns true.
    // Otherwise do nothing and return false.
    //
    bool splitTriangle(size_t tid) {
        for (int n = 0; n < 3; ++n) {
            int i = triangles[tid][n];
            int j = triangles[tid][(n + 1) % 3];
            int k = triangles[tid][(n + 2) % 3];

            // If edge [jk] intersects plane.
            int m = intersection(j, k);
            if (m != -1) {
                triangles.push_back({ i, j, m });
                triangles[tid] = { i, m, k };
                jud_plane.push_back(true);
                return true;
            }
        }
        return false;
    }
};

} // namespace vasco