#include <iostream>
#include "..\include\Element.h"

void Element::categorise()
{
    centroid = std::make_pair((vertices[0].first + vertices[3].first) / 2.0, (vertices[0].second + vertices[3].second) / 2.0);
    double di = trimmingCurve.projectionOfPoint(centroid);
    double r_in = sqrt(std::pow(vertices[3].first - centroid.first, 2));
    double r_out = sqrt(std::pow(vertices[3].first - centroid.first, 2) + std::pow(vertices[3].second - centroid.second, 2));
    computeTrimmedAndUntrimmedVertices();
    if (di < r_in) isTrimmed = true;
    else if (di > r_out) isTrimmed = false;
    else if (di > r_in && di < r_out)
    {
        if (trimmedVertices.size() == 0) isTrimmed = false;
        else isTrimmed = true;
    }
}

void Element::computeTrimmedAndUntrimmedVertices()
{
    if (trimmingCurve.isPointOutside(vertices[0])) trimmedVertices.push_back(vertices[0]);
    else untrimmedVertices.push_back(vertices[0]);
    if (trimmingCurve.isPointOutside(vertices[1])) trimmedVertices.push_back(vertices[1]);
    else untrimmedVertices.push_back(vertices[1]);
    if (trimmingCurve.isPointOutside(vertices[2])) trimmedVertices.push_back(vertices[2]);
    else untrimmedVertices.push_back(vertices[2]);
    if (trimmingCurve.isPointOutside(vertices[3])) trimmedVertices.push_back(vertices[3]);
    else untrimmedVertices.push_back(vertices[3]);
}

std::vector<std::vector<std::pair<double, double>>> Element::divideInTriangles()
{
    int amount = trimmedVertices.size();
    std::vector<std::vector<std::pair<double, double>>> triangles;
    switch (amount)
    {
    case 1:
        triangles = construct_3_triangles();
        break;
    case 2:
        triangles = construct_2_triangles();
        break;
    case 3:
        triangles = construct_1_triangle();
        break;
    default:
        break;
    }

    return triangles;
}

std::vector<std::vector<std::pair<double, double>>> Element::construct_3_triangles()
{
    double min_x = vertices[0].first;
    double max_x = vertices[3].first;
    double min_y = vertices[0].second;
    double max_y = vertices[3].second;
    std::pair<double, double> diagonal, vertex3, vertex4;

    bool XcoordinateOfFirstAppearsTwice = (untrimmedVertices[0].first == untrimmedVertices[1].first || untrimmedVertices[0].first == untrimmedVertices[2].first);
    bool YcoordinateOfFirstAppearsTwice = (untrimmedVertices[0].second == untrimmedVertices[1].second || untrimmedVertices[0].second == untrimmedVertices[2].second);
    bool XcoordinateOfSecondAppearsTwice = (untrimmedVertices[1].first == untrimmedVertices[0].first || untrimmedVertices[1].first == untrimmedVertices[2].first);
    bool YcoordinateOfSecondAppearsTwice = (untrimmedVertices[1].second == untrimmedVertices[0].second || untrimmedVertices[1].second == untrimmedVertices[2].second);
    if (XcoordinateOfFirstAppearsTwice && YcoordinateOfFirstAppearsTwice) diagonal = untrimmedVertices[0];
    else if (XcoordinateOfSecondAppearsTwice && YcoordinateOfSecondAppearsTwice) diagonal = untrimmedVertices[1];
    else diagonal = untrimmedVertices[2];

    for (int i = 0; i < 3; i++)
    {
        bool hasCommonYwithDiagonal = (untrimmedVertices[i].second == diagonal.second && untrimmedVertices[i] != diagonal);
        bool hasCommonXwithDiagonal = (untrimmedVertices[i].first == diagonal.first && untrimmedVertices[i] != diagonal);
        if (hasCommonYwithDiagonal) vertex3 = untrimmedVertices[i];
        if (hasCommonXwithDiagonal) vertex4 = untrimmedVertices[i];
    }

    std::pair<double, double> vertex1, vertex2;

    if (diagonal.first == min_x)
    {
        double s2 = max_x;
        double t = trimmingCurve.find_t_given_s(s2, min_y, max_y);
        vertex1 = std::make_pair(s2, t);

        if (diagonal.second == min_y)
        {
            double t1 = max_y;
            double s = trimmingCurve.find_s_given_t(t1, min_x, max_x);
            vertex2 = std::make_pair(s, t1);
        }
        else if (diagonal.second == max_y)
        {
            double t1 = min_y;
            double s = trimmingCurve.find_s_given_t(t1, min_x, max_x);
            vertex2 = std::make_pair(s, t1);
        }
    }
    else if (diagonal.first == max_x)
    {
        double s2 = min_x;
        double t = trimmingCurve.find_t_given_s(s2, min_y, max_y);
        vertex1 = std::make_pair(s2, t);

        if (diagonal.second == min_y)
        {
            double t1 = max_y;
            double s = trimmingCurve.find_s_given_t(t1, min_x, max_x);
            vertex2 = std::make_pair(s, t1);
        }
        else if (diagonal.second == max_y)
        {
            double t1 = min_y;
            double s = trimmingCurve.find_s_given_t(t1, min_x, max_x);
            vertex2 = std::make_pair(s, t1);
        }
    }

    std::vector<std::pair<double, double>> triangle1{vertex1, diagonal, vertex3};
    std::vector<std::pair<double, double>> triangle2{vertex1, diagonal, vertex2};
    std::vector<std::pair<double, double>> triangle3{vertex2, diagonal, vertex4};

    return std::vector<std::vector<std::pair<double, double>>>{triangle1, triangle2, triangle3};
}

std::vector<std::vector<std::pair<double, double>>> Element::construct_2_triangles()
{
    double min_x = vertices[0].first;
    double max_x = vertices[3].first;
    double min_y = vertices[0].second;
    double max_y = vertices[3].second;
    std::pair<double, double> vertex1, vertex2, vertex3;
    if (untrimmedVertices[0].first == untrimmedVertices[1].first)
    {
        if (untrimmedVertices[0].second > untrimmedVertices[1].second) vertex3 = untrimmedVertices[0];
        else vertex3 = untrimmedVertices[1];

        double t1 = max_y;
        double s = trimmingCurve.find_s_given_t(t1, min_x, max_x);
        vertex1 = std::make_pair(s, t1);

        t1 = min_y;
        s = trimmingCurve.find_s_given_t(t1, min_x, max_x);
        vertex2 = std::make_pair(s, t1);
    }
    else if (untrimmedVertices[0].second == untrimmedVertices[1].second)
    {
        if (untrimmedVertices[0].first > untrimmedVertices[1].first) vertex3 = untrimmedVertices[0];
        else vertex3 = untrimmedVertices[1];
        double s1 = max_x;
        double t = trimmingCurve.find_t_given_s(s1, min_y, max_y);
        vertex1 = std::make_pair(s1, t);

        s1 = min_x;
        t = trimmingCurve.find_t_given_s(s1, min_y, max_y);
        vertex2 = std::make_pair(s1, t);
    }

    std::vector<std::pair<double, double>> triangle1{vertex1, vertex3, vertex2};
    std::vector<std::pair<double, double>> triangle2{vertex2, untrimmedVertices[0], untrimmedVertices[1]};

    return std::vector<std::vector<std::pair<double, double>>>{triangle1, triangle2};
}

std::vector<std::vector<std::pair<double, double>>> Element::construct_1_triangle()
{
    double min_x = vertices[0].first;
    double max_x = vertices[3].first;
    double min_y = vertices[0].second;
    double max_y = vertices[3].second;

    double t1 = untrimmedVertices[0].second;
    double s = trimmingCurve.find_s_given_t(t1, min_x, max_x);
    std::pair<double, double> vertex1 = std::make_pair(s, t1);

    double s2 = untrimmedVertices[0].first;
    double t = trimmingCurve.find_t_given_s(s2, min_y, max_y);
    std::pair<double, double> vertex2 = std::make_pair(s2, t);

    std::vector<std::pair<double, double>> triangle1{vertex1, untrimmedVertices[0], vertex2};

    return std::vector<std::vector<std::pair<double, double>>>{triangle1};
}