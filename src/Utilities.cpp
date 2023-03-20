#include <iostream>
#include "..\include\Utilities.h"


namespace Utils
{

    template <class T>
    std::vector<T> createTensorProduct(std::vector<T> &vec1, std::vector<T> &vec2)
    {
        std::vector<T> tensor_product;
        for (int u = 0; u < vec1.size(); ++u)
        {
            for (int v = 0; v < vec2.size(); ++v)
            {
                tensor_product.push_back({vec1[u] * vec2[v]});
            }
        }

        return tensor_product;
    }

    template <class T>
    void mapValuesToDomain(std::vector<T> &GaussPoints, const T left, const T right)
    {
        assert(left <= right);
        for (int i = 0; i < GaussPoints.size(); ++i)
        {
            GaussPoints[i] = ((left * (1 - GaussPoints[i]) + right * (1 + GaussPoints[i])) / 2);
        }
    }

    template <class T>
    std::ostream &operator<<(std::ostream &os, const std::vector<T> &vec)
    {
        for (auto value : vec)
        {
            os << value << std::endl;
        }

        return os;
    }

    std::vector<std::pair<double, double>> GaussPointsAndWeightsQuad(int numberOfPoints, const double left, const double right)
    {
        assert(numberOfPoints > 0);
        assert(left <= right);

        std::vector<double> GaussPoints, GaussWeights;
        // Compute Gauss points in interval [0,1]
        switch (numberOfPoints)
        {
        case 1:
            GaussPoints = {0.0};
            GaussWeights = {2.0};
            break;
        case 2:
            GaussPoints = {-0.57735, 0.57735};
            GaussWeights = {1.0, 1.0};
            break;
        case 3:
            GaussPoints = {0.0, -0.774597, 0.774597};
            GaussWeights = {0.888889, 0.555556, 0.555556};
            break;
        case 4:
            GaussPoints = {-0.861136, -0.339981, 0.339981, 0.861136};
            GaussWeights = {0.347855, 0.652145, 0.652145, 0.347855};
            break;
        case 5:
            GaussPoints = {-0.90618, -0.538469, 0.0, 0.538469, 0.90618};
            GaussWeights = {0.236927, 0.478629, 0.568889, 0.478629, 0.236927};
            break;
        default:
            GaussPoints = {-0.90618, -0.538469, 0.0, 0.538469, 0.90618};
            GaussWeights = {0.236927, 0.478629, 0.568889, 0.478629, 0.236927};
            break;
        }
        mapValuesToDomain(GaussPoints, left, right);

        std::vector<std::pair<double, double>> GaussPointsAndWeights;
        for (int i = 0; i < GaussPoints.size(); ++i)
        {
            GaussPointsAndWeights.push_back(std::make_pair(GaussPoints[i], GaussWeights[i]));
        }

        return GaussPointsAndWeights;
    }

    std::vector<std::pair<double, double>> GaussPointsAndWeightsTriangle(double a, double b)
    {
        assert(a <= b);
        std::vector<double> GS_pts_temp{1.0 / 6.0, 2.0 / 3.0};
        std::vector<double> GS_wgts{1.0 / 3.0, 1.0 / 3.0};
        // Convert to interval [a,b]
        std::vector<double> GS_pts;
        for (int i = 0; i < GS_pts_temp.size(); ++i)
        {
            GS_pts.push_back(a * (1 - GS_pts_temp[i]) + b * GS_pts_temp[i]);
        }

        std::vector<std::pair<double, double>> GaussPointsAndWeights;
        for (int i = 0; i < GS_pts.size(); ++i)
        {
            GaussPointsAndWeights.push_back(std::make_pair(GS_pts[i], GS_wgts[i]));
        }

        return GaussPointsAndWeights;
    }

    template <typename T>
    T norm(const std::vector<T> &v)
    {
        T result = 0.0;
        for (size_t i = 0; i < v.size(); ++i)
        {
            result += v[i] * v[i];
        }

        return sqrt(result);
    }

} // end namespace Utilities