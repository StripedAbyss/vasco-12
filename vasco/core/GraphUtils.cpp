#include "GraphUtils.h"


namespace vasco {

std::vector<int> GeneralGraph_DArraySArraySpatVarying(
    int num_pixels,
    int num_labels,
    const std::vector<std::vector<int>>& data_value,
    const std::vector<std::vector<int>>& pixels_relations,
    const std::vector<int>& length_edges) noexcept
{
    std::vector<int> result(static_cast<size_t>(num_pixels), 0);

    // 数据项
    int* data = new (std::nothrow) int[static_cast<size_t>(num_pixels) * static_cast<size_t>(num_labels)];
    if (!data) {
        return result;
    }
    for (int i = 0; i < num_pixels; ++i) {
        for (int l = 0; l < num_labels; ++l) {
            data[i * num_labels + l] = data_value[static_cast<size_t>(i)][static_cast<size_t>(l)];
        }
    }

    // 光滑项（|l1-l2|==0 => 0, else 1）
    int* smooth = new (std::nothrow) int[static_cast<size_t>(num_labels) * static_cast<size_t>(num_labels)];
    if (!smooth) {
        delete[] data;
        return result;
    }
    for (int l1 = 0; l1 < num_labels; ++l1) {
        for (int l2 = 0; l2 < num_labels; ++l2) {
            smooth[l1 + l2 * num_labels] = (std::abs(l2 - l1) == 0) ? 0 : 1;
        }
    }

    try {
        int* label_cost = new int[static_cast<size_t>(num_labels)];
        for (int i = 0; i < num_labels; ++i) {
            label_cost[i] = 100000;
        }

        GCoptimizationGeneralGraph* gc = new GCoptimizationGeneralGraph(num_pixels, num_labels);
        gc->setDataCost(data);
        gc->setSmoothCost(smooth);
        gc->setLabelCost(label_cost);

        for (int i = 0; i < static_cast<int>(pixels_relations.size()); ++i) {
            gc->setNeighbors(pixels_relations[static_cast<size_t>(i)][0],
                             pixels_relations[static_cast<size_t>(i)][1],
                             length_edges[static_cast<size_t>(i)] / 10);
        }

        std::printf("\nBefore optimization energy is %lld", gc->compute_energy());
        gc->expansion(20);
        std::printf("\nAfter optimization energy is %lld", gc->compute_energy());

        std::cout << std::endl << "***********" << std::endl;
        std::cout << "DataEnergy: " << gc->giveDataEnergy() << std::endl;
        std::cout << "SmoothEnergy: " << gc->giveSmoothEnergy() << std::endl;
        std::cout << "LabelEnergy: " << gc->giveLabelEnergy() << std::endl;
        std::cout << "***********" << std::endl;

        for (int i = 0; i < num_pixels; ++i) {
            result[static_cast<size_t>(i)] = gc->whatLabel(i);
        }

        delete gc;
        delete[] label_cost;
    }
    catch (GCException& e) {
        e.Report();
    }

    delete[] smooth;
    delete[] data;

    return result;
}

std::vector<int> GeneralGraph_DArraySArraySpatVarying2(int num_pixels, int num_labels, const std::vector<std::vector<int>>& data_value,
    const std::vector<std::pair<int,int>>& pixels_relations,
    const std::vector<int>& length_edges) noexcept
{
    std::vector<int> result(static_cast<size_t>(num_pixels), 0);

    // 数据项
    int* data = new (std::nothrow) int[static_cast<size_t>(num_pixels) * static_cast<size_t>(num_labels)];
    if (!data) {
        return result;
    }
    for (int i = 0; i < num_pixels; ++i) {
        for (int l = 0; l < num_labels; ++l) {
            data[i * num_labels + l] = data_value[static_cast<size_t>(i)][static_cast<size_t>(l)];
        }
    }

    // 光滑项（|l1-l2|==0 => 0, else 1）
    int* smooth = new (std::nothrow) int[static_cast<size_t>(num_labels) * static_cast<size_t>(num_labels)];
    if (!smooth) {
        delete[] data;
        return result;
    }
    for (int l1 = 0; l1 < num_labels; ++l1) {
        for (int l2 = 0; l2 < num_labels; ++l2) {
            smooth[l1 + l2 * num_labels] = (std::abs(l2 - l1) == 0) ? 0 : 1;
        }
    }

    try {
        int* label_cost = new int[static_cast<size_t>(num_labels)];
        for (int i = 0; i < num_labels; ++i) {
            label_cost[i] = 100000;
        }

        GCoptimizationGeneralGraph* gc = new GCoptimizationGeneralGraph(num_pixels, num_labels);
        gc->setVerbosity(2);
        gc->setDataCost(data);
        gc->setSmoothCost(smooth);
        gc->setLabelCost(label_cost);

        /*for (int i = 0; i < static_cast<int>(pixels_relations.size()); ++i) {
            gc->setNeighbors(pixels_relations[static_cast<size_t>(i)].first,
                pixels_relations[static_cast<size_t>(i)].second,
                length_edges[static_cast<size_t>(i)]);
        }*/

		std::cout << "Before optimization energy is " << gc->compute_energy() << std::endl;
        std::cout << std::endl << "***********" << std::endl;
        std::cout << "DataEnergy: " << gc->giveDataEnergy() << std::endl;
        std::cout << "SmoothEnergy: " << gc->giveSmoothEnergy() << std::endl;
        std::cout << "LabelEnergy: " << gc->giveLabelEnergy() << std::endl;
        std::cout << "***********" << std::endl;

        gc->expansion(2);
		std::cout << "After optimization energy is " << gc->compute_energy() << std::endl;

        std::cout << std::endl << "***********" << std::endl;
        std::cout << "DataEnergy: " << gc->giveDataEnergy() << std::endl;
        std::cout << "SmoothEnergy: " << gc->giveSmoothEnergy() << std::endl;
        std::cout << "LabelEnergy: " << gc->giveLabelEnergy() << std::endl;
        std::cout << "***********" << std::endl;

        for (int i = 0; i < num_pixels; ++i) {
            result[static_cast<size_t>(i)] = gc->whatLabel(i);
        }

        delete gc;
        delete[] label_cost;
    }
    catch (GCException& e) {
        e.Report();
    }

    delete[] smooth;
    delete[] data;

    return result;
}

} // namespace vasco