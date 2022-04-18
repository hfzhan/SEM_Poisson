/*
 * generate 3D tetrahedral mesh file on $[0, a]^3$, 
 *     based on essential cube of spectral element method
 * 
 * argv[1]: size of domain $[0, a]^3$
 * argv[2]: number of equal partition in each direction
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <stdlib.h>
#include <iomanip>
#include <vector>

#define DIM 3

int main(int argc, char * argv[])
{
    std::ofstream outfile;
    outfile.open("cube.mesh");

    double a = atof(argv[1]);
    int n = atoi(argv[2]); // equal partition at each direction
    std::cout << "generate 3D tetrahedral mesh based on spectral element method\n" 
              << "truncated domain: [" << 0 << ", " << a << "]^3, each direction divided into " << n << " pieces " << std::endl;
    
    // std::cout << "n: " << n << std::endl;
    double h = a * 1. / n; // step size
    int n_point = (n+1) * (n+1) * (n+1); // number of points
    int n_edge_edge = n * (n+1) * (n+1); // number of edge locates on edge of cube, parallel to $x$ axis
    int n_edge_face = (n+1) * n * n; // number of edge locates on face of cube, orthogonal to $x$ axis
    int n_edge = n_edge_edge * 3 + n_edge_face * 3; // [parallel to one axis] + [orthogonal to one axis]
    int n_triangle_face = (n+1) * n * n * 2; // number of triangle locates on face of cube, orthogonal to $x$ axis
    int n_triangle_inner = n * n * n * 4;// number of triangle locates in cube
    int n_triangle = n_triangle_face * 3 + n_triangle_inner; // [faces of cube] + [inner triangle in cube]
    int n_tetrahedron = n * n * n * 5; // each cube has 5 tetrahedron
    int flag = 0; // flag of boundary element

    // arrange point in row-major order, which means consider 3D multiindex (i1, i2, i3) corresponds to point number (i1 * n + i2) * n + i3
    // while multiindex (i1, i2, i3) corresponds to point (-a, -a, -a) * (i1, i2, i3) * h
    // output point info
    outfile << n_point << '\n'; // number of points
    double h_tmp = 1.0 / n;
    for (int i1 = 0; i1 < n+1; ++i1)
        for (int i2 = 0; i2 < n+1; ++i2)
            for (int i3 = 0; i3 < n+1; ++i3)
                outfile << i1*h << '\t' << i2*h << '\t' << i3*h << '\n';
                // outfile << i1*h_tmp << '\t' << i2*h_tmp << '\t' << i3*h_tmp << '\n';
    std::cerr << "have output point coordinates\n";

    // output 0-dimensional geometry
    // attach boundary label 1 to the point locating on the faces of $[-a, a]^3$
    outfile << '\n' << n_point << '\n';
    // int point[n_point][3]; // point with multiindex (i1, i2, i3) = (i1, i2, i3) * h + (1, 1, 1) * (-a)
    std::vector<std::vector<int> > point;
    point.resize(n_point);
    for (int i = 0; i < n_point; ++i)
        point[i].resize(3);
    // std::cerr << "initialize point\n";
    int flag_point[n_point]; // flag of point
    for (int i3 = 0; i3 < n+1; ++i3)
        for (int i2 = 0; i2 < n+1; ++i2)
            for (int i1 = 0; i1 < n+1; ++i1){
                int index_point = (i1*(n+1) + i2)*(n+1) + i3;
                flag = 0; // flag of inner point
                if (i1 == 0 || i1 == n || i2 == 0 || i2 == n || i3 == 0 || i3 == n)
                    flag = 1; // flag of boundary point
                point[index_point][0] = i1;
                point[index_point][1] = i2;
                point[index_point][2] = i3;
                flag_point[index_point] = flag;
                // outfile << index_point << "\n\t"
                //         << 1 << '\t' << index_point << "\n\t"
                //         << 1 << '\t' << index_point << "\n\t"
                //         << flag << '\n';
            }
    for (int index_point = 0; index_point < n_point; ++index_point)
        outfile << index_point << "\n\t"
                << 1 << '\t' << index_point << "\n\t"
                << 1 << '\t' << index_point << "\n\t"
                << flag_point[index_point] << '\n';
    std::cerr << "have output 0-dimensional info\n";

    // divide edges into 6 classes, classes 0 to 2: edges of cube, parallel to x, y, z axis, on face of cube, classes 3 to 5: orthogonal to x, y, z axis
    // denote the endpoints of edge as (x1, y1, z1) and (x2, y2, z2)
    // represent each edge by multiindex (i1, i2, i3), corresponds to the multiindex of the point (min(x1,x2), min(y1,y2), min(z1,z2))
    // output 1-dimension geometry
    outfile << '\n' << n_edge << '\n'; // number of edges
    // output 1-dimensional geometry locates on edges of cube
    // int endpoint[n_edge][2]; // endpoint of each edge
    std::vector<std::vector<int > > endpoint;
    endpoint.resize(n_edge);
    for (int i = 0; i < n_edge; ++i)
        endpoint[i].resize(2);
    // std::cerr << "initialize endpoint\n";
    // int flag_edge[n_edge];
    std::vector<int> flag_edge;
    flag_edge.resize(n_edge);
    // std::cerr << "initialize flag_edge\n";
    for (int i1 = 0; i1 < n; ++i1)
        for (int i2 = 0; i2 < n+1; ++i2)
            for (int i3 = 0; i3 < n+1; ++i3){
                // these three class of edge share the same boundary flag for the same multiindex (i1, i2, i3)
                // std::cout << "i1 = " << i1 << ", i2 = " << i2 << ", i3 = " << i3 << '\n';
                flag = 0; // flag of inner edge
                if (i2 == 0 || i2 == n || i3 == 0 || i3 == n)
                    flag = 1; // flag of boundary edge
                // output edges in class 1 with multiindex (i1, i2, i3): locates on edge of cube, parallel to x axis
                int index_edge = (i1*(n+1) + i2)*(n+1) + i3;
                int index_point_start = (i1*(n+1) + i2)*(n+1) + i3;
                int index_point_end = ((i1+1)*(n+1) + i2)*(n+1) + i3;
                if (index_edge >= n_edge_edge)
                    std::cout << "error! i1 = " << i1 << ", i2 = " << i2 << ", i3 = " << i3 << ", index_edge = " << index_edge << '\n';
                if ((i1 + i2 + i3) % 2 == 0){ // exchange two endpoints
                    int tmp = index_point_start;
                    index_point_start = index_point_end;
                    index_point_end = tmp;
                }
                endpoint[index_edge][0] = index_point_start;
                endpoint[index_edge][1] = index_point_end;
                flag_edge[index_edge] = flag;
                // outfile << index_edge << "\n\t"
                //         << 2 << '\t' << index_point_start << '\t' << index_point_end << "\n\t"
                //         << 2 << '\t' << index_point_start << '\t' << index_point_end << "\n\t"
                //         << flag << '\n';
                // output edges in class 2 with multiindex (i2, i1, i3): locates on edge of cube, parallel to y axis
                index_edge = (i2*n + i1)*(n+1) + i3 + n_edge_edge;
                index_point_start = (i2*(n+1) + i1)*(n+1) + i3;
                index_point_end = (i2*(n+1) + i1+1)*(n+1) + i3;
                if (index_edge >= n_edge_edge*2)
                    std::cout << "error! i1 = " << i2 << ", i2 = " << i1 << ", i3 = " << i3 << ", index_edge = " << index_edge << '\n';
                if ((i1 + i2 + i3) % 2 == 0){ // exchange two endpoints
                    int tmp = index_point_start;
                    index_point_start = index_point_end;
                    index_point_end = tmp;
                }
                endpoint[index_edge][0] = index_point_start;
                endpoint[index_edge][1] = index_point_end;
                flag_edge[index_edge] = flag;
                // outfile << index_edge << "\n\t"
                //         << 2 << '\t' << index_point_start << '\t' << index_point_end << "\n\t"
                //         << 2 << '\t' << index_point_start << '\t' << index_point_end << "\n\t"
                //         << flag << '\n';
                // output edges in class 3 with multiindex (i3, i2, i1): locates on edge of cube, parallel to z axis
                index_edge = (i3*(n+1) + i2)*n + i1 + n_edge_edge*2;
                index_point_start = (i3*(n+1) + i2)*(n+1) + i1;
                index_point_end = (i3*(n+1) + i2)*(n+1) + i1+1;
                if (index_edge >= n_edge_edge*3)
                    std::cout << "error! i1 = " << i3 << ", i2 = " << i2 << ", i3 = " << i1 << ", index_edge = " << index_edge << '\n';
                if ((i1 + i2 + i3) % 2 == 0){ // exchange two endpoints
                    int tmp = index_point_start;
                    index_point_start = index_point_end;
                    index_point_end = tmp;
                }
                endpoint[index_edge][0] = index_point_start;
                endpoint[index_edge][1] = index_point_end;
                flag_edge[index_edge] = flag;
                // outfile << index_edge << "\n\t"
                //         << 2 << '\t' << index_point_start << '\t' << index_point_end << "\n\t"
                //         << 2 << '\t' << index_point_start << '\t' << index_point_end << "\n\t"
                //         << flag << '\n';
            }
    // std::cerr << "assign edge locates on edges of cube\n";
    // output 1-dimensional geometry locates on faces of cubes
    for (int i1 = 0; i1 < n+1; ++i1)
        for (int i2 = 0; i2 < n; ++i2)
            for (int i3 = 0; i3 < n; ++i3){
                // these three class of edge share the same boundary flag for the same multiindex (i1, i2, i3)
                flag = 0; // flag of inner edge
                if (i1 == 0 || i1 == n)
                    flag = 1; // flag of boundary edge
                // output edges in class 4 with multiindex (i1, i2, i3): locates on edge of cube, orthogonal to x axis
                int index_point_start, index_point_end;
                int index_edge = (i1*n + i2)*n + i3 + n_edge_edge*3;
                if ((i1 + i2 + i3) % 2 == 0){ // edge tpye 1: point (i1, i2, i3) is an endpoint
                    index_point_start = (i1*(n+1) + i2)*(n+1) + i3; // point with smalled index
                    index_point_end = (i1*(n+1) + i2+1)*(n+1) + i3+1;// point with larger index
                    if ((i1 + i2) % 2 == 1){
                        int tmp = index_point_start;
                        index_point_start = index_point_end;
                        index_point_end = tmp;
                    }
                }
                else{ // edge type 2: point (i1, i2, i3) is not an endpoint
                    index_point_start = (i1*(n+1) + i2+1)*(n+1) + i3; // point with smaller z coordinate
                    index_point_end = (i1*(n+1) + i2)*(n+1) + i3+1; // point with larger z coordinate
                    if ((i1 + i2) % 2 == 0){
                        int tmp = index_point_start;
                        index_point_start = index_point_end;
                        index_point_end = tmp;
                    }
                }
                endpoint[index_edge][0] = index_point_start;
                endpoint[index_edge][1] = index_point_end;
                flag_edge[index_edge] = flag;
                // outfile << index_edge << "\n\t"
                //         << 2 << '\t' << index_point_start << '\t' << index_point_end << "\n\t"
                //         << 2 << '\t' << index_point_start << '\t' << index_point_end << "\n\t"
                //         << flag << '\n';
                // output edges in class 5 with multiindex (i2, i1, i3): locates on edge of cube, orthogonal to y axis
                index_edge = (i2*(n+1) + i1)*n + i3 + n_edge_edge*3 + n_edge_face;
                if ((i1 + i2 + i3) % 2 == 0){ // edge tpye 1: point (i2, i1, i3) is an endpoint
                    index_point_start = (i2*(n+1) + i1)*(n+1) + i3; // point with smalled index
                    index_point_end = ((i2+1)*(n+1) + i1)*(n+1) + i3+1;// point with larger index
                    if ((i2 + i1) % 2 == 1){
                        int tmp = index_point_start;
                        index_point_start = index_point_end;
                        index_point_end = tmp;
                    }
                }
                else{ // edge type 2: point (i2, i1, i3) is not an endpoint
                    index_point_start = (i2*(n+1) + i1)*(n+1) + i3+1; // point with larger z coordinate
                    index_point_end = ((i2+1)*(n+1) + i1)*(n+1) + i3; // point with smaller z coordinate
                    if ((i2 + i1) % 2 == 1){
                        int tmp = index_point_start;
                        index_point_start = index_point_end;
                        index_point_end = tmp;
                    }
                }
                endpoint[index_edge][0] = index_point_start;
                endpoint[index_edge][1] = index_point_end;
                flag_edge[index_edge] = flag;
                // outfile << index_edge << "\n\t"
                //         << 2 << '\t' << index_point_start << '\t' << index_point_end << "\n\t"
                //         << 2 << '\t' << index_point_start << '\t' << index_point_end << "\n\t"
                //         << flag << '\n';
                // output edges in class 6 with multiindex (i3, i2, i1): locates on edge of cube, orthogonal to z axis
                index_edge = (i3*n + i2)*(n+1) + i1 + n_edge_edge*3 + n_edge_face*2;
                if ((i1 + i2 + i3) % 2 == 0){ // edge tpye 1: point (i3, i2, i1) is an endpoint
                    index_point_start = ((i3+1)*(n+1) + i2+1)*(n+1) + i1; // point with smalled index
                    index_point_end = (i3*(n+1) + i2)*(n+1) + i1;// point with larger index
                    if ((i3 + i1) % 2 == 1){
                        int tmp = index_point_start;
                        index_point_start = index_point_end;
                        index_point_end = tmp;
                    }
                }
                else{ // edge type 2: point (i2, i1, i3) is not an endpoint
                    index_point_start = ((i3+1)*(n+1) + i2)*(n+1) + i1; // point with smaller y coordinate
                    index_point_end = (i3*(n+1) + i2+1)*(n+1) + i1; // point with larger y coordinate
                    if ((i3 + i1) % 2 == 1){
                        int tmp = index_point_start;
                        index_point_start = index_point_end;
                        index_point_end = tmp;
                    }
                }
                endpoint[index_edge][0] = index_point_start;
                endpoint[index_edge][1] = index_point_end;
                flag_edge[index_edge] = flag;
                // outfile << index_edge << "\n\t"
                //         << 2 << '\t' << index_point_start << '\t' << index_point_end << "\n\t"
                //         << 2 << '\t' << index_point_start << '\t' << index_point_end << "\n\t"
                //         << flag << '\n';
            }
    int n_edge_fixedstart[n_point]; // [i]: number of edges whose start point is i-th point
    for (int index_point = 0; index_point < n_point; ++index_point)
        n_edge_fixedstart[index_point] = 0;
    // int index_edge_fixedstart[n_point][18]; // [i][]: the index of edge whose start point is i-th point
    std::vector<std::vector<int> > index_edge_fixedstart;
    index_edge_fixedstart.resize(n_point);
    for (int i = 0; i < n_point; ++i)
        index_edge_fixedstart[i].resize(18);
    for (int index_edge = 0; index_edge < n_edge; ++index_edge){
        outfile << index_edge << "\n\t"
                << 2 << '\t' << endpoint[index_edge][0] << '\t' << endpoint[index_edge][1] << "\n\t"
                << 2 << '\t' << endpoint[index_edge][0] << '\t' << endpoint[index_edge][1] << "\n\t"
                << flag_edge[index_edge] << '\n';
        index_edge_fixedstart[endpoint[index_edge][0]][n_edge_fixedstart[endpoint[index_edge][0]]] = index_edge;
        n_edge_fixedstart[endpoint[index_edge][0]]++;
    }
    std::cerr << "have output 1-dimensional info\n";
    // std::cout << "test for index_edge_fixedstart\n";
    // for (int i = 0; i < n_point; ++i){
    //     std::cout << "point index: " << i << " with " << n_edge_fixedstart[i] << " edges\n";
    //     for (int ind = 0; ind < n_edge_fixedstart[i]; ++ind){
    //         int ind_index = index_edge_fixedstart[i][ind];
    //         std::cout << "edge index: " << ind_index << " with endpoint " << endpoint[ind_index][0] << ", " << endpoint[ind_index][1] << '\n';
    //     }
    // }

    // output 2-dimensional geometry
    // output 2-dimensional geometry locates on face of cubes
    // int point_triangle[n_triangle][3]; // point indices of triangles with correspondence: [0]: (0, 0); [1]: (1, 0); [2]: (0, 1)
    // int edge_triangle[n_triangle][3]; // edge indices of triangles with correspondence: [0]: (0, 0); [1]: (1, 0); [2]: (0, 1)
    // int flag_triangle[n_triangle]; // flag of triangle
    std::vector<std::vector<int> > point_triangle;
    std::vector<std::vector<int> > edge_triangle;
    std::vector<int> flag_triangle;
    point_triangle.resize(n_triangle);
    edge_triangle.resize(n_triangle);
    flag_triangle.resize(n_triangle);
    for (int i = 0; i < n_triangle; ++i){
        point_triangle[i].resize(3);
        edge_triangle[i].resize(3);
    }
    // traverse face of cubes, assign triangle info for clase 1 to 3
    for (int i1 = 0; i1 < n+1; ++i1)
        for (int i2 = 0; i2 < n; ++i2)
            for (int i3 = 0; i3 < n; ++i3){
                // output face is class 1 with multiindex (i1, i2, i3): locates on face of cubes, orthogonal to x axis
                int index_triangle[2]; // [0]: lower triangle; [1]: upper triangle
                int index_hypertenuse;
                int index_point_original[2]; // [0]: lower; [1]: upper 
                index_triangle[0] = ((i1*n + i2)*n + i3) * 2;
                index_triangle[1] = ((i1*n + i2)*n + i3) * 2 + 1; 
                index_hypertenuse = (i1*n + i2)*n + i3 + n_edge_edge*3;
                flag = 0;
                if (i1 == 0 || i1 == n)
                    flag = 1;
                // find the indices of original points on lower/upper reference triangle
                if ((i1 + i2 + i3) % 2 == 0){ // (i1, i2, i3) is an endpoint of hypertenuse, then two original points are (i1, i2+1, i3), (i1, i2, i3+1)
                    index_point_original[0] = (i1*(n+1) + i2+1)*(n+1) + i3;
                    index_point_original[1] = (i1*(n+1) + i2)*(n+1) + i3+1;
                }
                else{ // (i1, i2, i3) is an original point in reference triangle, then (i1, i2+1, i3+1) is also an original point
                    index_point_original[0] = (i1*(n+1) + i2)*(n+1) + i3;
                    index_point_original[1] = (i1*(n+1) + i2+1)*(n+1) + i3+1;
                }
                if (i3 % 2 == 1){ // if i3 is odd, exchange two original points, and corresponding edges
                    int tmp = index_point_original[0];
                    index_point_original[0] = index_point_original[1];
                    index_point_original[1] = tmp;
                }
                // find the indices of edges on lower/upper reference triangle
                int index_edge_local[2][2]; // [0][0]: lower x; [0][1]: lower y; [1][0]: upper x; [1][1]: upper y
                for (int ind_lu = 0; ind_lu < 2; ++ind_lu)
                    for (int ind_xy = 0; ind_xy < 2; ++ind_xy)
                        // for (int index_edge = 0; index_edge < n_edge; ++index_edge)
                        //     if (endpoint[index_edge][0] == index_point_original[ind_lu] && // go from original point
                        //         endpoint[index_edge][1] == endpoint[index_hypertenuse][ind_xy]){ // x: left endpoint of hypertenuse; y: right endpoint of hypertenuse
                        //         index_edge_local[ind_lu][ind_xy] = index_edge;
                        //         break;
                        //     }
                        for (int ind = 0; ind < n_edge_fixedstart[index_point_original[ind_lu]]; ++ind){
                            int ind_edge = index_edge_fixedstart[index_point_original[ind_lu]][ind]; // index of edge whose start point is l/u original point
                            if (endpoint[ind_edge][1] == endpoint[index_hypertenuse][ind_xy]){ // its end point is the same as corresponding point of hypertenuse
                                index_edge_local[ind_lu][ind_xy] = ind_edge;
                            //     std::cout << "original point: " << index_point_original[ind_lu] << ", end point: " << endpoint[index_hypertenuse][ind_xy] << '\n';
                            //     std::cout << "index of hypertenuse: " << index_hypertenuse << ", start point: " << endpoint[index_hypertenuse][0] << ", end point: " << endpoint[index_hypertenuse][1] << '\n';
                            //     std::cout << "edge index: " << ind_edge << ", start point: " << endpoint[ind_edge][0] << ", end point: " << endpoint[ind_edge][1] << '\n';
                            }
                        }
                // assign triangle info for class 1
                for (int ind_lu = 0; ind_lu < 2; ++ind_lu){
                    point_triangle[index_triangle[ind_lu]][0] = index_point_original[ind_lu];
                    edge_triangle[index_triangle[ind_lu]][0] = index_hypertenuse;
                    for (int ind_xy = 1; ind_xy <= 2; ++ind_xy){
                        point_triangle[index_triangle[ind_lu]][ind_xy] = endpoint[index_hypertenuse][ind_xy-1];
                        edge_triangle[index_triangle[ind_lu]][ind_xy] = index_edge_local[ind_lu][2-ind_xy];
                    }
                    flag_triangle[index_triangle[ind_lu]] = flag;
                }
                // exchange index such that traverse each info anti-clockwise
                // int ind_exchange;
                // if ((i2 + i3) % 2 == 0)
                //     ind_exchange = 1;
                // else
                //     ind_exchange = 0;
                // int tmp_exchange = point_triangle[index_triangle[ind_exchange]][1];              
                // point_triangle[index_triangle[ind_exchange]][1] = point_triangle[index_triangle[ind_exchange]][2];
                // point_triangle[index_triangle[ind_exchange]][2] = tmp_exchange;
                // tmp_exchange = edge_triangle[index_triangle[ind_exchange]][1];              
                // edge_triangle[index_triangle[ind_exchange]][1] = edge_triangle[index_triangle[ind_exchange]][2];
                // edge_triangle[index_triangle[ind_exchange]][2] = tmp_exchange;
                // find face info class 2 with multiindex (i2, i1, i3): locates on face of cubes, orthogonal to y axis
                index_triangle[0] = ((i2*(n+1) + i1)*n + i3) * 2 + n_triangle_face;
                index_triangle[1] = ((i2*(n+1) + i1)*n + i3) * 2 + 1 + n_triangle_face; 
                index_hypertenuse = (i2*(n+1) + i1)*n + i3 + n_edge_edge*3 + n_edge_face;
                // find the indices of original points on lower/upper reference triangle
                if ((i1 + i2 + i3) % 2 == 0){ // (i2, i1, i3) is an endpoint of hypertenuse, then two original points are (i2+1, i2, i3), (i2, i1, i3+1)
                    index_point_original[0] = ((i2+1)*(n+1) + i1)*(n+1) + i3;
                    index_point_original[1] = (i2*(n+1) + i1)*(n+1) + i3+1;
                }
                else{ // (i2, i1, i3) is an original point in reference triangle, then (i2+1, i1, i3+1) is also an original point
                    index_point_original[0] = (i2*(n+1) + i1)*(n+1) + i3;
                    index_point_original[1] = ((i2+1)*(n+1) + i1)*(n+1) + i3+1;
                }
                if (i3 % 2 == 1){ // if i3 is odd, exchange two original points, and corresponding edges
                    int tmp = index_point_original[0];
                    index_point_original[0] = index_point_original[1];
                    index_point_original[1] = tmp;
                }
                // find the indices of edges on lower/upper reference triangle
                for (int ind_lu = 0; ind_lu < 2; ++ind_lu)
                    for (int ind_xy = 0; ind_xy < 2; ++ind_xy)
                        // for (int index_edge = 0; index_edge < n_edge; ++index_edge)
                        //     if (endpoint[index_edge][0] == index_point_original[ind_lu] && // go from original point
                        //         endpoint[index_edge][1] == endpoint[index_hypertenuse][ind_xy]){ // x: start point of hypertenuse; y: end point of hypertenuse
                        //         index_edge_local[ind_lu][ind_xy] = index_edge;
                        //         break;
                        //     }
                        for (int ind = 0; ind < n_edge_fixedstart[index_point_original[ind_lu]]; ++ind)
                            if (endpoint[index_edge_fixedstart[index_point_original[ind_lu]][ind]][1] == endpoint[index_hypertenuse][ind_xy])
                                index_edge_local[ind_lu][ind_xy] = index_edge_fixedstart[index_point_original[ind_lu]][ind];
                // assign triangle info for class 2
                for (int ind_lu = 0; ind_lu < 2; ++ind_lu){
                    point_triangle[index_triangle[ind_lu]][0] = index_point_original[ind_lu];
                    edge_triangle[index_triangle[ind_lu]][0] = index_hypertenuse;
                    for (int ind_xy = 0; ind_xy < 2; ++ind_xy){
                        point_triangle[index_triangle[ind_lu]][ind_xy+1] = endpoint[index_hypertenuse][ind_xy];
                        edge_triangle[index_triangle[ind_lu]][ind_xy+1] = index_edge_local[ind_lu][1-ind_xy];
                    }
                    flag_triangle[index_triangle[ind_lu]] = flag;
                }
                // find face info class 3 with multiindex (i3, i2, i1): locates on face of cubes, orthogonal to z axis
                index_triangle[0] = ((i3*n + i2)*(n+1) + i1) * 2 + n_triangle_face*2;
                index_triangle[1] = ((i3*n + i2)*(n+1) + i1) * 2 + 1 + n_triangle_face*2; 
                index_hypertenuse = (i3*n + i2)*(n+1) + i1 + n_edge_edge*3 + n_edge_face*2;
                // find the indices of original points on lower/upper reference triangle
                if ((i1 + i2 + i3) % 2 == 0){ // (i3, i2, i1) is an endpoint of hypertenuse, then two original points are (i3, i2+1, i1), (i3+1, i2, i1)
                    index_point_original[0] = (i3*(n+1) + i2+1)*(n+1) + i1;
                    index_point_original[1] = ((i3+1)*(n+1) + i2)*(n+1) + i1;
                }
                else{ // (i3, i2, i1) is an original point in reference triangle, then (i3+1, i2+1, i1) is also an original point
                    index_point_original[0] = (i3*(n+1) + i2)*(n+1) + i1;
                    index_point_original[1] = ((i3+1)*(n+1) + i2+1)*(n+1) + i1;
                }
                if (i3 % 2 == 1){ // if i3 is odd, exchange two original points, and corresponding edges
                    int tmp = index_point_original[0];
                    index_point_original[0] = index_point_original[1];
                    index_point_original[1] = tmp;
                }
                // find the indices of edges on lower/upper reference triangle
                for (int ind_lu = 0; ind_lu < 2; ++ind_lu)
                    for (int ind_xy = 0; ind_xy < 2; ++ind_xy)
                        // for (int index_edge = 0; index_edge < n_edge; ++index_edge)
                        //     if (endpoint[index_edge][0] == index_point_original[ind_lu] && // go from original point
                        //         endpoint[index_edge][1] == endpoint[index_hypertenuse][ind_xy]){ // x: start point of hypertenuse; y: end point of hypertenuse
                        //         index_edge_local[ind_lu][ind_xy] = index_edge;
                        //         break;
                        //     }
                        for (int ind = 0; ind < n_edge_fixedstart[index_point_original[ind_lu]]; ++ind)
                            if (endpoint[index_edge_fixedstart[index_point_original[ind_lu]][ind]][1] == endpoint[index_hypertenuse][ind_xy])
                                index_edge_local[ind_lu][ind_xy] = index_edge_fixedstart[index_point_original[ind_lu]][ind];
                // assign triangle info for class 3
                for (int ind_lu = 0; ind_lu < 2; ++ind_lu){
                    point_triangle[index_triangle[ind_lu]][0] = index_point_original[ind_lu];
                    edge_triangle[index_triangle[ind_lu]][0] = index_hypertenuse;
                    for (int ind_xy = 0; ind_xy < 2; ++ind_xy){
                        point_triangle[index_triangle[ind_lu]][ind_xy+1] = endpoint[index_hypertenuse][ind_xy];
                        edge_triangle[index_triangle[ind_lu]][ind_xy+1] = index_edge_local[ind_lu][1-ind_xy];
                    }
                    flag_triangle[index_triangle[ind_lu]] = flag;
                }
            }
    // traverse cube, assign triangle info for triangle in cubes
    for (int i1 = 0; i1 < n; ++i1)
        for (int i2 = 0; i2 < n; ++i2)
            for (int i3 = 0; i3 < n; ++i3){ // the point (i1, i2, i3), node of cube with least coordinates
                int index_triangle_start = ((i1*n + i2)*n + i3) * 4 + n_triangle_face*3;
                int index_hypertenuse[6];
                index_hypertenuse[0] = (i1*n + i2)*n + i3 + n_edge_edge*3; // hypertenuse locates on faces orthogonal to x axis
                index_hypertenuse[1] = ((i1+1)*n + i2)*n + i3 + n_edge_edge*3; // [0]: x -> y; [1]: O -> z
                index_hypertenuse[2] = (i1*(n+1) + i2)*n + i3 + n_edge_edge*3 + n_edge_face; // hypertenuse locates on faces orthogonal to y axis
                index_hypertenuse[3] = (i1*(n+1) + i2+1)*n + i3 + n_edge_edge*3 + n_edge_face; // [2]: x -> z; [3]: O -> y
                index_hypertenuse[4] = (i1*n + i2)*(n+1) + i3 + n_edge_edge*3 + n_edge_face*2; // hypertenuse locates on faces orthogonal to z axis
                index_hypertenuse[5] = (i1*n + i2)*(n+1) + i3+1 + n_edge_edge*3 + n_edge_face*2; // [4]: O -> x; [5]: y -> z
                // for (int ind = 0; ind < 6; ++ind)
                //     outfile << "index_hypertenuse[" << ind << "] = " << index_hypertenuse[ind] << '\n';
                if (i1 % 2 == 1){
                    int tmp = index_hypertenuse[0];
                    index_hypertenuse[0] = index_hypertenuse[1];
                    index_hypertenuse[1] = tmp;
                }
                if (i2 % 2 == 1){
                    int tmp = index_hypertenuse[2];
                    index_hypertenuse[2] = index_hypertenuse[3];
                    index_hypertenuse[3] = tmp;
                }
                if (i3 % 2 == 1){
                    int tmp = index_hypertenuse[4];
                    index_hypertenuse[4] = index_hypertenuse[5];
                    index_hypertenuse[5] = tmp;
                }
                int index_point_local[4] = {endpoint[index_hypertenuse[4]][0], endpoint[index_hypertenuse[4]][1],
                                            endpoint[index_hypertenuse[5]][0], endpoint[index_hypertenuse[5]][1]}; // [0]: O; [1]: x; [2]: y; [3]: z
                // for (int ind = 0; ind < 4; ++ind)
                //     outfile << "index_point_local[" << ind << "] = " << index_point_local[ind] << '\n';
                // // assign triangle info for clase 4
                // triangle Oxy
                int index_triangle_now;
                index_triangle_now = index_triangle_start;
                point_triangle[index_triangle_now][0] = index_point_local[0];
                point_triangle[index_triangle_now][1] = index_point_local[1];
                point_triangle[index_triangle_now][2] = index_point_local[2];
                edge_triangle[index_triangle_now][0] = index_hypertenuse[0];
                edge_triangle[index_triangle_now][1] = index_hypertenuse[3];
                edge_triangle[index_triangle_now][2] = index_hypertenuse[4];
                flag_triangle[index_triangle_now] = 0;
                // triangle Oxz
                index_triangle_now = index_triangle_start+1;
                point_triangle[index_triangle_now][0] = index_point_local[0];
                point_triangle[index_triangle_now][1] = index_point_local[1];
                point_triangle[index_triangle_now][2] = index_point_local[3];
                edge_triangle[index_triangle_now][0] = index_hypertenuse[2];
                edge_triangle[index_triangle_now][1] = index_hypertenuse[1];
                edge_triangle[index_triangle_now][2] = index_hypertenuse[4];
                flag_triangle[index_triangle_now] = 0;
                // triangle xyz
                index_triangle_now = index_triangle_start+2;
                point_triangle[index_triangle_now][0] = index_point_local[1];
                point_triangle[index_triangle_now][1] = index_point_local[2];
                point_triangle[index_triangle_now][2] = index_point_local[3];
                edge_triangle[index_triangle_now][0] = index_hypertenuse[5];
                edge_triangle[index_triangle_now][1] = index_hypertenuse[2];
                edge_triangle[index_triangle_now][2] = index_hypertenuse[0];
                flag_triangle[index_triangle_now] = 0;
                // triangle Oyz
                index_triangle_now = index_triangle_start+3;
                point_triangle[index_triangle_now][0] = index_point_local[0];
                point_triangle[index_triangle_now][1] = index_point_local[2];
                point_triangle[index_triangle_now][2] = index_point_local[3];
                edge_triangle[index_triangle_now][0] = index_hypertenuse[5];
                edge_triangle[index_triangle_now][1] = index_hypertenuse[1];
                edge_triangle[index_triangle_now][2] = index_hypertenuse[3];
                flag_triangle[index_triangle_now] = 0;
            }
    // output triangle info
    outfile << '\n' << n_triangle << '\n'; // number of triangle
    for (int index_triangle = 0; index_triangle < n_triangle; ++index_triangle)
        outfile << index_triangle << "\n\t"
                << 3 << '\t' << point_triangle[index_triangle][0] << '\t' << point_triangle[index_triangle][1] << '\t' << point_triangle[index_triangle][2] << "\n\t"
                << 3 << '\t' << edge_triangle[index_triangle][0] << '\t' << edge_triangle[index_triangle][1] << '\t' << edge_triangle[index_triangle][2] << "\n\t"
                << flag_triangle[index_triangle] << "\n";
    std::cerr << "have output 2-dimensional info\n";

    // output 3-dimensional info
    // int point_tetrahedron[n_tetrahedron][4];
    // int face_tetrahedron[n_tetrahedron][4];
    // int flag_tetrahedron[n_tetrahedron];
    std::vector<std::vector<int> > point_tetrahedron;
    std::vector<std::vector<int> > face_tetrahedron;
    std::vector<int> flag_tetrahedron;
    // std::cerr << "name quantities\n";
    point_tetrahedron.resize(n_tetrahedron);
    face_tetrahedron.resize(n_tetrahedron);
    flag_tetrahedron.resize(n_tetrahedron);
    // std::cerr << "resize quantities\n";
    for (int i = 0; i < n_tetrahedron; ++i){
        point_tetrahedron[i].resize(4);
        face_tetrahedron[i].resize(4);
    }
    // std::cerr << "initialize tetrahedron info\n";
    // traverse cube, assign tetrahedron info
    // for (int sum = 0; sum <= (n-1)*3; ++sum) // traverse cube according to summation of index
    //     for (int i1 = 0; i1 <= std::min(n-1, sum); ++i1)
    //         for (int i2 = 0; i2 <= std::min(n-1, sum-i1); ++i2){
                // int i3 = sum - i1 -i2;
                // if (i3 >= n) continue;
    for (int i1 = 0; i1 < n; ++i1)
        for (int i2 = 0; i2 < n; ++i2)
            for (int i3 = 0; i3 < n; ++i3){
                int index_tetrahedron_start = ((i1*n + i2)*n + i3) * 5;
                // find 0-dimensional info
                int index_point_local[8]; // nodes of cube. on reference cube, [0]: (0, 0, 0); [1]: (1, 0, 0); [2]: (0, 1, 0); [4]: (0, 0, 1)
                int d1 = i1 % 2, d2 = i2 % 2, d3 = i3 % 2;
                for (int j1 = 0; j1 < 2; ++j1) // consider mirror tranformation
                    for (int j2 = 0; j2 < 2; ++j2)
                        for (int j3 = 0; j3 < 2; ++j3){
                            // index_point_local[j1 + 2*j2 + 4*j3] = ((i1+((j1+d1)%2))*(n+1) + i2+((j2+d2)%2))*(n+1) + i3+((j3+d3)%2);
                            int shift1 = (i1 + j1) % 2, shift2 = (i2 + j2) % 2, shift3 = (i3 + j3) % 2;
                            index_point_local[(shift1*2 + shift2)*2 + shift3] = ((i1+j1)*(n+1) + i2+j2)*(n+1) + i3+j3;
                        }
                // find 2-dimensional info
                int index_triangle_local[16]; // triangle of cube, take mirror transformation into consideration
                for (int j = 0; j < 2; ++j) // triangle locates on faces of cube, orthogonal to x, y, z axis respectively, [even]: lower; [odd]: upper
                    for (int d = 0; d < 2; ++d){
                        // int i1_m = (i1 + d1) % 2, i2_m = (i2 + d2) % 2, i3_m = (i3 + d3) % 2; // m for modified, considering the mirror transfromation
                        // int i1_am = (i1 + d + d1) % 2, i2_am = (i2 + d + d2) % 2, i3_am = (i3 + d + d3) % 2; // a for addtion
                        // index_triangle_local[0 + d*2 + j] = ((((i1+d+d1)%2)*n + i2)*n + i3) * 2 + j; // least coordinate: (i1+d, i2, i3), orthogonal to x axis
                        // index_triangle_local[4 + d*2 + j] = ((i1*(n+1) + (i2+d+d2%2))*n + i3_m) * 2 + n_triangle_face + j; // (i1, i2+d, i3), y axis
                        // index_triangle_local[8 + d*2 + j] = ((i1*n + i2)*(n+1) + i3_am) * 2 + n_triangle_face*2 + j; // (i1, i2, i3+d), z axis
                        int shift1 = (i1 + d) % 2, shift2 = (i2 + d) % 2, shift3 = (i3 + d) % 2;
                        index_triangle_local[0 + shift1*2 + j] = (((i1+d)*n + i2)*n + i3) * 2 + j;
                        index_triangle_local[4 + shift2*2 + j] = ((i1*(n+1) + i2+d)*n + i3) * 2 + j + n_triangle_face;
                        index_triangle_local[8 + shift3*2 + j] = ((i1*n + i2)*(n+1) + i3+d) * 2 + j + n_triangle_face*2;
                }
                for (int ind = 0; ind < 4; ++ind) // triangle consisting of hypertenuses
                    index_triangle_local[12 + ind] = ((i1*n + i2)*n + i3) * 4 + n_triangle_face*3 + ind;
                // assign tetrahedron
                int index_tetrahedron_now;
                // set [0] to be the tetrahedron consisting of hypertenuses
                index_tetrahedron_now = index_tetrahedron_start;
                point_tetrahedron[index_tetrahedron_now][0] = index_point_local[6];
                point_tetrahedron[index_tetrahedron_now][1] = index_point_local[0];
                point_tetrahedron[index_tetrahedron_now][2] = index_point_local[3];
                point_tetrahedron[index_tetrahedron_now][3] = index_point_local[5];
                face_tetrahedron[index_tetrahedron_now][0] = index_triangle_local[14];
                face_tetrahedron[index_tetrahedron_now][1] = index_triangle_local[15];
                face_tetrahedron[index_tetrahedron_now][2] = index_triangle_local[13];
                face_tetrahedron[index_tetrahedron_now][3] = index_triangle_local[12];
                // traverse tetrahedron anti-clockwise
                index_tetrahedron_now = index_tetrahedron_start+1;
                point_tetrahedron[index_tetrahedron_now][0] = index_point_local[1];
                point_tetrahedron[index_tetrahedron_now][1] = index_point_local[0];
                point_tetrahedron[index_tetrahedron_now][2] = index_point_local[3];
                point_tetrahedron[index_tetrahedron_now][3] = index_point_local[5];
                face_tetrahedron[index_tetrahedron_now][0] = index_triangle_local[14];
                face_tetrahedron[index_tetrahedron_now][1] = index_triangle_local[10];
                face_tetrahedron[index_tetrahedron_now][2] = index_triangle_local[5];
                face_tetrahedron[index_tetrahedron_now][3] = index_triangle_local[1];
                index_tetrahedron_now = index_tetrahedron_start+2;
                point_tetrahedron[index_tetrahedron_now][0] = index_point_local[4];
                point_tetrahedron[index_tetrahedron_now][1] = index_point_local[6];
                point_tetrahedron[index_tetrahedron_now][2] = index_point_local[0];
                point_tetrahedron[index_tetrahedron_now][3] = index_point_local[5];
                face_tetrahedron[index_tetrahedron_now][0] = index_triangle_local[13];
                face_tetrahedron[index_tetrahedron_now][1] = index_triangle_local[4];
                face_tetrahedron[index_tetrahedron_now][2] = index_triangle_local[2];
                face_tetrahedron[index_tetrahedron_now][3] = index_triangle_local[9];
                index_tetrahedron_now = index_tetrahedron_start+3;
                point_tetrahedron[index_tetrahedron_now][0] = index_point_local[7];
                point_tetrahedron[index_tetrahedron_now][1] = index_point_local[6];
                point_tetrahedron[index_tetrahedron_now][2] = index_point_local[3];
                point_tetrahedron[index_tetrahedron_now][3] = index_point_local[5];
                face_tetrahedron[index_tetrahedron_now][0] = index_triangle_local[15];
                face_tetrahedron[index_tetrahedron_now][1] = index_triangle_local[11];
                face_tetrahedron[index_tetrahedron_now][2] = index_triangle_local[3];
                face_tetrahedron[index_tetrahedron_now][3] = index_triangle_local[7];
                index_tetrahedron_now = index_tetrahedron_start+4;
                point_tetrahedron[index_tetrahedron_now][0] = index_point_local[2];
                point_tetrahedron[index_tetrahedron_now][1] = index_point_local[6];
                point_tetrahedron[index_tetrahedron_now][2] = index_point_local[0];
                point_tetrahedron[index_tetrahedron_now][3] = index_point_local[3];
                face_tetrahedron[index_tetrahedron_now][0] = index_triangle_local[12];
                face_tetrahedron[index_tetrahedron_now][1] = index_triangle_local[0];
                face_tetrahedron[index_tetrahedron_now][2] = index_triangle_local[6];
                face_tetrahedron[index_tetrahedron_now][3] = index_triangle_local[8];
                // assign flag of tetrahedron
                for (int ind = 0; ind < 5; ++ind){ // set boundary element if one of the faces of tetrahedron locates on boundary
                    flag = 0;
                    for (int ind_face = 0; ind_face < 4; ++ind_face)
                        if (flag_triangle[face_tetrahedron[index_tetrahedron_start+ind][ind_face]])
                            flag = 1;
                    flag_tetrahedron[index_tetrahedron_start+ind] = flag;
                }
            }
    // output 3-dimensional info
    int order[2][4]; // [0]: original; [1]: exchange 1, 2 and 3
    order[0][0] = 0; order[0][1] = 1; order[0][2] = 2; order[0][3] = 3;
    order[1][0] = 0; order[1][1] = 2; order[1][2] = 1; order[1][3] = 3;
    outfile << '\n' << n_tetrahedron << '\n'; //number of tetrahedron
    for (int index_tetrahedron = 0; index_tetrahedron < n_tetrahedron; ++index_tetrahedron){
        double p[4][3]; // [i][j]: i-th point in this tetrahedron, j-th coordinate
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 3; ++j)
                p[i][j] = point[point_tetrahedron[index_tetrahedron][i]][j] * h;
        for (int i = 1; i < 4; ++i) // calculate the difference
            for (int j = 0; j < 3; ++j)
                p[i][j] -= p[0][j];
        double volume;
        volume = p[1][0]*p[2][1]*p[3][2] + p[1][1]*p[2][2]*p[3][0] + p[1][2]*p[2][0]*p[3][1] - p[1][0]*p[2][2]*p[3][1]
            - p[1][1]*p[2][0]*p[3][2] - p[1][2]*p[2][1]*p[3][0];
        int ind_order = 0;
        if (volume < 0) ind_order = 1;
        ind_order = 0;
        outfile << index_tetrahedron << "\n\t4";
        for (int ind = 0; ind < 4; ++ind)
            outfile << '\t' << point_tetrahedron[index_tetrahedron][order[ind_order][ind]];
        outfile << "\n\t4";
        for (int ind = 0; ind < 4; ++ind)
            outfile << '\t' << face_tetrahedron[index_tetrahedron][order[ind_order][ind]];
        outfile << "\n\t" << flag_tetrahedron[index_tetrahedron] << '\n';
    }
    std::cerr << "have output 3-dimensional info\n";
    
    outfile.close();
}
