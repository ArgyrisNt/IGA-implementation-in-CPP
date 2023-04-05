#include <iostream>
#include "..\include\BsplineEntity.h"

void BsplineEntity::knotInsertion(KnotVector<double> &vector, std::vector<Vertex<double>> &points, const double newKnot)
{
    KnotVector<double> new_knotvector = vector;
    int span = vector.findSpanOfValue(newKnot);
    new_knotvector.insert(span + 1, newKnot);

    std::vector<Vertex<double>> newPoints;
    std::vector<double> newWeights;
    for (int j = span - vector.getDegree() + 1; j <= span; ++j)
    {
        double alpha = (newKnot - vector(j)) / (vector(j + vector.getDegree()) - vector(j));
        double w1 = vector.getWeights()[j - 1];
        double w2 = vector.getWeights()[j];
        double coord_x = (1.0 - alpha) * points[j - 1].x * w1 + alpha * points[j].x * w2;
        double coord_y = (1.0 - alpha) * points[j - 1].y * w1 + alpha * points[j].y * w2;
        double new_weight = (1.0 - alpha) * w1 + alpha * w2;
        newWeights.push_back(new_weight);
        newPoints.push_back({coord_x / new_weight, coord_y / new_weight});
    }

    std::vector<Vertex<double>> FinalNewPoints;
    std::vector<double> FinalNewWeights;
    for (int i = 0; i <= span - vector.getDegree(); ++i)
    {
        FinalNewPoints.push_back(points[i]);
        FinalNewWeights.push_back(vector.getWeights()[i]);
    }
    for (int ii = 0; ii < newPoints.size(); ++ii)
    {
        FinalNewPoints.push_back(newPoints[ii]);
        FinalNewWeights.push_back(newWeights[ii]);
    }
    for (int ii = span; ii < points.size(); ++ii)
    {
        FinalNewPoints.push_back(points[ii]);
        FinalNewWeights.push_back(vector.getWeights()[ii]);
    }
    
    new_knotvector.setWeights(FinalNewWeights);
    points = FinalNewPoints;
    vector = new_knotvector;
}

void BsplineEntity::plotControlPoints(const std::string &filename)
{
    std::ofstream plotCtrlPts(filename);
    plotCtrlPts << "X,Y\n";
    for (auto it = controlPoints.begin(); it != controlPoints.end(); ++it)
    {
        plotCtrlPts << (*it).x << "," << (*it).y << "\n";
    }
    plotCtrlPts.close();
}