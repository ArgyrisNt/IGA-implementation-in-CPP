#include <iostream>
#include "..\include\BsplineEntity.h"

void BsplineEntity::knotInsertion(KnotVector<double> &vector, std::vector<Vertex<double>> &points, const double newKnot)
{
    KnotVector<double> new_knotvector = vector;
    int span = vector.findSpanOfValue(newKnot);
    new_knotvector.insert(span + 1, newKnot);

    std::vector<Vertex<double>> newPoints;
    for (int j = span - vector.getDegree() + 1; j <= span; ++j)
    {
        double alpha = (newKnot - vector(j)) / (vector(j + vector.getDegree()) - vector(j));
        double coord_x = (1.0 - alpha) * points[j - 1].x + alpha * points[j].x;
        double coord_y = (1.0 - alpha) * points[j - 1].y + alpha * points[j].y;
        newPoints.push_back({coord_x, coord_y});
    }

    std::vector<Vertex<double>> FinalNewPoints;
    for (int i = 0; i <= span - vector.getDegree(); ++i)
    {
        FinalNewPoints.push_back(points[i]);
    }
    for (int ii = 0; ii < newPoints.size(); ++ii)
    {
        FinalNewPoints.push_back(newPoints[ii]);
    }
    for (int ii = span; ii < points.size(); ++ii)
    {
        FinalNewPoints.push_back(points[ii]);
    }

    points = FinalNewPoints;
    vector = new_knotvector;
}