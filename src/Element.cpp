#include <iostream>
#include "..\include\Element.h"

void Element::categorise()
{
    centroid = Vertex<double>((vertices[0].x + vertices[3].x) / 2.0, (vertices[0].y + vertices[3].y) / 2.0);
    double di = trimmingCurve.projectionOfPoint(centroid);
    double r_in = sqrt(std::pow(vertices[3].x - centroid.x, 2));
    double r_out = sqrt(std::pow(vertices[3].x - centroid.x, 2) + std::pow(vertices[3].y - centroid.y, 2));
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

std::vector<Triangle<double>> Element::divideInTriangles()
{
    int amount = trimmedVertices.size();
    std::vector<Triangle<double>> triangles;
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

std::vector<Triangle<double>> Element::construct_3_triangles()
{
    double min_x = vertices[0].x;
    double max_x = vertices[3].x;
    double min_y = vertices[0].y;
    double max_y = vertices[3].y;
    Vertex<double> diagonal, vertex3, vertex4;

    bool XcoordOfFirstAppearsTwice = (almostEqual(untrimmedVertices[0].x, untrimmedVertices[1].x) || almostEqual(untrimmedVertices[0].x, untrimmedVertices[2].x));
    bool YcoordOfFirstAppearsTwice = (almostEqual(untrimmedVertices[0].y, untrimmedVertices[1].y) || almostEqual(untrimmedVertices[0].y, untrimmedVertices[2].y));
    bool XcoordOfSecondAppearsTwice = (almostEqual(untrimmedVertices[1].x, untrimmedVertices[0].x) || almostEqual(untrimmedVertices[1].x, untrimmedVertices[2].x));
    bool YcoordOfSecondAppearsTwice = (almostEqual(untrimmedVertices[1].y, untrimmedVertices[0].y) || almostEqual(untrimmedVertices[1].y, untrimmedVertices[2].y));
    if (XcoordOfFirstAppearsTwice && YcoordOfFirstAppearsTwice) diagonal = untrimmedVertices[0];
    else if (XcoordOfSecondAppearsTwice && YcoordOfSecondAppearsTwice) diagonal = untrimmedVertices[1];
    else diagonal = untrimmedVertices[2];

    for (int i = 0; i < 3; i++)
    {
        bool hasCommonYwithDiagonal = (almostEqual(untrimmedVertices[i].y, diagonal.y) && untrimmedVertices[i] != diagonal);
        bool hasCommonXwithDiagonal = (almostEqual(untrimmedVertices[i].x, diagonal.x) && untrimmedVertices[i] != diagonal);
        if (hasCommonYwithDiagonal) vertex3 = untrimmedVertices[i];
        if (hasCommonXwithDiagonal) vertex4 = untrimmedVertices[i];
    }

    Vertex<double> vertex1, vertex2;

    if (almostEqual(diagonal.x, min_x))
    {
        double s2 = max_x;
        double t = trimmingCurve.find_t_given_s(s2, min_y, max_y);
        vertex1.set(s2, t);

        if (almostEqual(diagonal.y, min_y))
        {
            double t1 = max_y;
            double s = trimmingCurve.find_s_given_t(t1, min_x, max_x);
            vertex2.set(s, t1);
        }
        else if (almostEqual(diagonal.y, max_y))
        {
            double t1 = min_y;
            double s = trimmingCurve.find_s_given_t(t1, min_x, max_x);
            vertex2.set(s, t1);
        }
    }
    else if (almostEqual(diagonal.x, max_x))
    {
        double s2 = min_x;
        double t = trimmingCurve.find_t_given_s(s2, min_y, max_y);
        vertex1.set(s2, t);

        if (almostEqual(diagonal.y, min_y))
        {
            double t1 = max_y;
            double s = trimmingCurve.find_s_given_t(t1, min_x, max_x);
            vertex2.set(s, t1);
        }
        else if (almostEqual(diagonal.y, max_y))
        {
            double t1 = min_y;
            double s = trimmingCurve.find_s_given_t(t1, min_x, max_x);
            vertex2.set(s, t1);
        }
    }

    Triangle<double> triangle1{vertex1, diagonal, vertex3};
    Triangle<double> triangle2{vertex1, diagonal, vertex2};
    Triangle<double> triangle3{vertex2, diagonal, vertex4};

    return std::vector<Triangle<double>>{triangle1, triangle2, triangle3};
}

std::vector<Triangle<double>> Element::construct_2_triangles()
{
    double min_x = vertices[0].x;
    double max_x = vertices[3].x;
    double min_y = vertices[0].y;
    double max_y = vertices[3].y;
    Vertex<double> vertex1, vertex2, vertex3;
    if (almostEqual(untrimmedVertices[0].x, untrimmedVertices[1].x))
    {
        if (untrimmedVertices[0].y > untrimmedVertices[1].y) vertex3 = untrimmedVertices[0];
        else vertex3 = untrimmedVertices[1];

        double t1 = max_y;
        double s = trimmingCurve.find_s_given_t(t1, min_x, max_x);
        vertex1.set(s, t1);

        t1 = min_y;
        s = trimmingCurve.find_s_given_t(t1, min_x, max_x);
        vertex2.set(s, t1);
    }
    else if (almostEqual(untrimmedVertices[0].y, untrimmedVertices[1].y))
    {
        if (untrimmedVertices[0].x > untrimmedVertices[1].x) vertex3 = untrimmedVertices[0];
        else vertex3 = untrimmedVertices[1];
        double s1 = max_x;
        double t = trimmingCurve.find_t_given_s(s1, min_y, max_y);
        vertex1.set(s1, t);

        s1 = min_x;
        t = trimmingCurve.find_t_given_s(s1, min_y, max_y);
        vertex2.set(s1, t);
    }

    Triangle<double> triangle1(vertex1, vertex3, vertex2);
    Triangle<double> triangle2(vertex2, untrimmedVertices[0], untrimmedVertices[1]);

    return std::vector<Triangle<double>>{triangle1, triangle2};
}

std::vector<Triangle<double>> Element::construct_1_triangle()
{
    double min_x = vertices[0].x;
    double max_x = vertices[3].x;
    double min_y = vertices[0].y;
    double max_y = vertices[3].y;

    double t1 = untrimmedVertices[0].y;
    double s = trimmingCurve.find_s_given_t(t1, min_x, max_x);
    Vertex<double> vertex1(s, t1);

    double s2 = untrimmedVertices[0].x;
    double t = trimmingCurve.find_t_given_s(s2, min_y, max_y);
    Vertex<double> vertex2(s2, t);

    Triangle<double> triangle1(vertex1, untrimmedVertices[0], vertex2);

    return std::vector<Triangle<double>>{triangle1};
}