/*
 * visualization of sem solution
 *
 * argv[1]: file name of mesh
 * argv[2]: file name of sem solution
 * argv[3]: maximum number of local refinement
 * argv[4]: tolerance for break criterion of local refinement
 */

#include <iostream>
#include <fstream>
#include <iomanip>

#include <omp.h>

#include <AFEPack/Geometry.h>
#include <AFEPack/TemplateElement.h>
#include <AFEPack/FEMSpace.h>
#include <AFEPack/HGeometry.h>

#define PI (4.0*atan(1.0))
#define DIM 3

const int N_GlobalRefine = 2;

std::vector<int> n_transform_local;
std::vector<std::vector<int> > transform_local;
std::vector<std::vector<double> > weight_transform_local;
int M;

double calc_volume_tetrahedron(const double * v0, const double * v1, const double * v2, const double * v3)
{
	return fabs(((v1[0] - v0[0]) * (v2[1] - v0[1]) * (v3[2] - v0[2])
                 + (v1[1] - v0[1]) * (v2[2] - v0[2]) * (v3[0] - v0[0])
                 + (v1[2] - v0[2]) * (v2[0] - v0[0]) * (v3[1] - v0[1])
                 - (v1[0] - v0[0]) * (v2[2] - v0[2]) * (v3[1] - v0[1])
                 - (v1[1] - v0[1]) * (v2[0] - v0[0]) * (v3[2] - v0[2])
                 - (v1[2] - v0[2]) * (v2[1] - v0[1]) * (v3[0] - v0[0]))/6.);
}

double calc_volume_triangle(const double * v0, const double * v1, const double * v2)
{
    double p1[DIM], p2[DIM];
    for (int i = 0; i < DIM; ++i){
        p1[i] = v1[i] - v0[i];
        p2[i] = v2[i] - v0[i];
    }
    return .5 * fabs(sqrt(pow(p1[1] * p2[2] - p1[2] * p2[1], 2) +
                          pow(p1[0] * p2[2] - p1[2] * p2[0], 2) +
                          pow(p1[0] * p2[1] - p1[1] * p2[0], 2)));
}

double calc_err(const std::vector<AFEPack::Point<DIM> >& vertex, const std::vector<AFEPack::Point<DIM> >& vertex_p,
                const std::vector<AFEPack::Point<DIM> >& QPoint, const std::vector<double>& Weight,
                const std::vector<double>& sol_local, const std::vector<double>& sol_sem_local);

int calc_binomial_coefficient(int n, int m)
{// calculate binomial coefficient $\binom{n}{m}$
    int val = 1;
    if (n > 0 && n >= m){
        if (n-m < m) m = n-m;
        for (int i = 0; i < m; ++i){val *= (n-i); val /= (i+1);}}
    return val;
}

struct Multiindex
{// DIM dimensional multiindex
    int index[DIM];
    bool operator <= (const Multiindex& multiindex){ // overload relation <= for multiindex
        for (int i = 0; i < DIM; ++i) if (this->index[i] > multiindex.index[i]) return false;
        return true;
    }
    bool operator == (const Multiindex& multiindex){ // overload relation <= for multiindex
        for (int i = 0; i < DIM; ++i) if (this->index[i] != multiindex.index[i]) return false;
        return true;
    }
    bool operator != (const Multiindex& multiindex){ // overload relation != for multiindex
        for (int i = 0; i < DIM; ++i) if (this->index[i] != multiindex.index[i]) return true;
        return false;}
    void operator = (const Multiindex& multiindex){ // overload operation = for multiindex
        for (int i = 0; i < DIM; ++i) this->index[i] = multiindex.index[i];}
    Multiindex operator + (const Multiindex& multiindex){ // overload operation + for multiindex
        Multiindex ans;
        for (int i = 0; i < DIM; ++i) ans.index[i] = this->index[i] + multiindex.index[i];
        return ans;}
    Multiindex operator - (const Multiindex& multiindex){ // overload operation + for multiindex
        Multiindex ans;
        for (int i = 0; i < DIM; ++i) ans.index[i] = this->index[i] - multiindex.index[i];
        return ans;}
    Multiindex operator * (int n){ // overload operation * for multiindex
        Multiindex ans;
        for (int i = 0; i < DIM; ++i) ans.index[i] = this->index[i] * n;
        return ans;}
    int sum(){ // summation of all components
        int ans;
        for (int i = 0; i < DIM; ++i) ans += this->index[i];
        return ans;}
    int n_nonzero(){ // number of nonzero components
        int ans = 0;
        for (int i = 0; i < DIM; ++i) if (this->index[i] != 0) ans++;
        return ans;}
};

static Multiindex Unitary_Multiindex[DIM];

int calc_binomial_coefficient_multiindex(Multiindex alpha, Multiindex beta)
{// calculate binomial coefficient for multiindex $\binom{\alpha}{\beta}$
    int val = 1;
    for (int ind = 0; ind < DIM; ++ind) // for $\beta$ or $\alpha-\beta$ with negative component, return 0
        if (beta.index[ind] < 0 || beta.index[ind] > alpha.index[ind])
            return 0;
    for (int ind = 0; ind < DIM; ++ind)
        val *= calc_binomial_coefficient(alpha.index[ind], beta.index[ind]);
    return val;
}

class Correspondence
{/* correpondence between
  *   multiindex and number of multiindex in lexicographical order
  */
    int M; // truncation order
    int number_index; // number of multiindex $|\alpha|\leq M$
    Multiindex *number_to_index; // number_to_index[i] = i-th multiindex in lexicographical order
    int **number_multiindex; /* number_multiindex[i=0:DIM-1]: number of multiindices with $(i+1)$ components and summation less than some number
                              *     number_multiindex[i=0:DIM-1][j=0:M]: number of multiindices $\beta=(\beta_k)_{k=0}^i\leq j$
                              */
public:
    // Multiindex unitary_multiindex[DIM];
    void init(int order_truncate); // initialize with truncation order
    int index2number(Multiindex multiindex); // return the number of given multiindex in lexicographical order
    Multiindex number2index(int number_index){ // return the (number_index)-th multiindex
        return number_to_index[number_index];}
    int n_index(){ // return number of multiindex $|\alpha|\leq M$, where M is truncation order
        return number_multiindex[DIM-1][M];}
    int n_index_begin(int order){ // return begin number of multiindex for $|\alpha|=order$ in lexicographical order
        if (order == 0) return 0;
        else return number_multiindex[DIM-1][order-1];}
    int n_index_end(int order){ // return begin number of multiindex for $|\alpha|=order+1$ in lexicographical order, which is also the number of multiindex $|\alpha|\leq order$
        return number_multiindex[DIM-1][order];}
} correspondence;

void Correspondence::init(int order_truncate)
{// initialize correspondence with given truncation order
    // unitary multiindex
    for (int i = 0; i < DIM; ++i){
        for (int j = 0; j < DIM; ++j)
            Unitary_Multiindex[i].index[j] = 0;
        Unitary_Multiindex[i].index[i] = 1;
    }
    M = order_truncate;
    // calculate number_multiindex
    number_multiindex = new int* [DIM]; // number_multiindex[i=0:DIM-1]: number of multiindices with $(i+1)$ components and summation less than some number
    for (int i = 0; i < DIM; ++i)
        number_multiindex[i] = new int [M+1];
    for (int i = 0; i < DIM; ++i) // number_multiindex[i=0:DIM-1][j=0:M]: number of multiindices $\beta=(\beta_k)_{k=0}^i\leq j$
        for (int j = 0; j <= M; ++j){
            number_multiindex[i][j] = calc_binomial_coefficient(i+j, i);
            if (j > 0) number_multiindex[i][j] += number_multiindex[i][j-1];
        }
    number_index = number_multiindex[DIM-1][M];
    // calculate i-th multiindex
    number_to_index = new Multiindex [number_index];
    for (int i = 0; i < number_index; ++i){
        int tmp_i = i+1;
        int sum = M;
        for (int j = 0; j < DIM; ++j){ // traverse summation of multiindex $|(\alpha_k)_{k=j}^{DIM-1}|$
            int k = sum - 1;
            for (; k >= 0; k -= 1)
                if (number_multiindex[DIM-1-j][k] < tmp_i)
                    break;
            if (j > 0) number_to_index[i].index[j-1] = sum - (k+1); // now step = 1
            sum = k + 1;
            if (k >= 0) tmp_i -= number_multiindex[DIM-1-j][k];
        }
        number_to_index[i].index[DIM-1] = sum;
    }
}

int Correspondence::index2number(Multiindex multiindex)
{// calculate the number of given multiindex
    int num = 0; // default for zero multiindex
    for (int i = 0; i < DIM; ++i){
        int sum = 0; // summation of last several components of multiindex $|(\alpha_k)_{k=i}^{DIM-1}|$
        for (int j = i; j < DIM; ++j)
            sum += multiindex.index[j];
        if (sum == 0) break;
        else num += number_multiindex[DIM-1-i][sum-1];
    }
    return num;
}

double calc_generalized_jacobi_polynomial(int alpha, int beta, int k, double x)
{// calculate k-th generalized jocobi polynomial with index (alpha, beta) at point x, where $alpha, beta >= -1$
    if (k == 0) return 1.0;
    if (alpha + beta == -2)
        if (k == 1) return x;
        else return 0.25 * (x-1) * (x+1) * calc_generalized_jacobi_polynomial(1, 1, k-2, x);
    if (alpha == -1)
        return (k+beta) * (x-1) / (k*2) * calc_generalized_jacobi_polynomial(1, beta, k-1, x);
    if (beta == -1)
        return (k+alpha) * (x+1) / (k*2) * calc_generalized_jacobi_polynomial(alpha, 1, k-1, x);
    double ans = 0.0, tmp_power = 1.0;
    for (int j = 0; j <= k; ++j){
        double factor = 1.0;
        for (int i = 0; i < k-j; ++i)
            factor *= (alpha+j+1 + i) / (i+1);
        for (int i = 0; i < j; ++i)
            factor *= (k+alpha+beta+1 + i) / (i+1);
        ans += factor * tmp_power;
        tmp_power *= (x - 1) * 0.5;
    }
    return ans;
}

double calc_coefficient_a(int alpha, int beta, int ind, int k)
{// calculate coefficient $a_{ind, k}^{alpha, beta}$
    if (k < 0) return 0;
    switch (ind){
    case 1:
        if (alpha == -1 && beta == -1)
            switch (k){
            case 0: return 1;
            case 1: return 4;
            case 2: return 0.5;
            }
        if (k == 0 && alpha + beta != -2)
            return 2.0 / (alpha + beta + 2);
        return 2.0 * (k+1) * (k+alpha+beta+1) / ((2*k+alpha+beta+1) * (2*k+alpha+beta+2));
    case 2:
        if (0 <= k && k <= 2 && alpha == -1 && beta == -1) return 0;
        if (k == 0 && alpha + beta != -2)
            return (beta - alpha) / (alpha + beta + 2.0);
        return (pow(beta,2)-pow(alpha,2)) / ((2*k+alpha+beta) * (2*k+alpha+beta+2.0));
    case 3:
        if (alpha == -1 && beta == -1){
            if (k == 0) return 0;
            if (k == 1) return 1;
            if (k == 2) return 0;
        }
        if (k == 0 && alpha+beta != -2)
            return 0;
        return 2.0 * (k+alpha) * (k+beta) / ((2*k+alpha+beta) * (2*k+alpha+beta+1));
    }
    return 0;
}

int main(int argc, char * argv[])
{
    // read mesh from argv[1]
    HGeometryTree<DIM> h_tree;
    h_tree.readMesh(argv[1]);
    IrregularMesh<DIM> *irregular_mesh;
    irregular_mesh = new IrregularMesh<DIM>;
    irregular_mesh->reinit(h_tree);
    irregular_mesh->semiregularize();
    irregular_mesh->regularize(false);
    RegularMesh<DIM> &mesh = irregular_mesh->regularMesh();
    int n_geometry[DIM+1];
    for (int ind = 0; ind <= DIM; ++ind)
        n_geometry[ind] = mesh.n_geometry(ind);
    // calculate refencen point for each geometry
    std::vector<std::vector<AFEPack::Point<DIM> > > point_ref_mesh; // [i = 0:3], barycenter of i dimensional geometry
    point_ref_mesh.resize(4);
    for (int ind = 0; ind <= DIM; ++ind){
        point_ref_mesh[ind].resize(mesh.n_geometry(ind));
        for (int i = 0; i < mesh.n_geometry(ind); ++i){
            AFEPack::Point<DIM> point_tmp = mesh.point(mesh.geometry(ind, i).vertex(0));
            for (int indt = 1; indt <= ind; ++indt)
                point_tmp += mesh.point(mesh.geometry(ind, i).vertex(indt));
            point_tmp /= (ind + 1);
            point_ref_mesh[ind][i] = point_tmp;
        }
    }

    
    // generate finite element space
    TemplateGeometry<DIM> template_geometry;
    CoordTransform<DIM, DIM> coord_transform;
    TemplateDOF<DIM> template_dof;
    BasisFunctionAdmin<double, DIM, DIM> basis_function;
    template_geometry.readData("tetrahedron.tmp_geo");
    coord_transform.readData("tetrahedron.crd_trs");
    template_dof.reinit(template_geometry);
    template_dof.readData("tetrahedron.2.tmp_dof");
    basis_function.reinit(template_dof);
    basis_function.readData("tetrahedron.2.bas_fun");
    std::vector<TemplateElement<double, DIM, DIM> > template_element(1);
    template_element[0].reinit(template_geometry,
                               template_dof,
                               coord_transform,
                               basis_function);
    FEMSpace<double, DIM> fem_space(mesh, template_element);
    int n_element = mesh.n_geometry(DIM);
    fem_space.element().resize(n_element);
    for (int i = 0; i < n_element; ++i)
        fem_space.element(i).reinit(fem_space, i, 0);
    fem_space.buildElement();
    fem_space.buildDof();
    fem_space.buildDofBoundaryMark();
    FEMSpace<double, DIM>::ElementIterator the_element = fem_space.beginElement();
    FEMSpace<double, DIM>::ElementIterator end_element = fem_space.endElement();
    
    
    // read SEM solution for argv[2]
    std::ifstream inputfile;
    int n_dof_total;
    inputfile.open(argv[2]);
    inputfile >> M;
    inputfile >> n_dof_total;
    Vector<double> sol_sem(n_dof_total);
    for (int i = 0; i < n_dof_total; ++i)
        inputfile >> sol_sem(i);
    inputfile.close();

    
    // initialize
    int n_dof = fem_space.n_dof();
    int n_fem_element = fem_space.n_element();
    correspondence.init(M);
    int n_index = correspondence.n_index();
    std::cout << "Initialize correspondence of multiindex with truncation order " << M << ", find " << n_index << " multiindices\n";
    std::vector<int> n_dof_geometry(DIM+1); // number of dof location on each dimensional geoemtry
    for (int ind = 0; ind <= DIM; ++ind) // number of dof on ind dimensional geometry
        n_dof_geometry[ind] = calc_binomial_coefficient(M-1, ind);
    for (int ind = 0; ind <= DIM; ++ind)
        std::cout << "n_dof_geometry[" << ind << "] = " << n_dof_geometry[ind] << ",\t";
    std::cout << '\n';
    std::cout << "the total dof is " << n_dof_total << '\n';
    
    
    // read quadrature info
    int n_quad_accuracy[3] = {38, 20, 20}; // [i]: accuracy of i dimensional quadrature rule
    int n_q_point[3] = {20, 88, 448}; // [i]: number of quadrature point in i dimensional rule
    std::vector<AFEPack::Point<DIM> > QPoint; // 3-d quadrature point
    std::vector<std::vector<std::vector<double> > > QPoint_Barycenter; // 1 & 2-d quadrature point in barycenter coordinate
    QPoint.resize(n_q_point[DIM-1]);
    QPoint_Barycenter.resize(DIM-1);
    for (int ind = 0; ind < DIM-1; ++ind){
        QPoint_Barycenter[ind].resize(n_q_point[ind]);
        for (int i = 0; i < n_q_point[ind]; ++i)
            QPoint_Barycenter[ind][i].resize(ind+2);
    }
    std::vector<std::vector<double> > Weight;
    Weight.resize(DIM);
    for (int ind = 0; ind < DIM; ++ind)
        Weight[ind].resize(n_q_point[ind]);
    for (int ind = 0; ind < DIM; ++ind){
        switch (ind){
        case 0: inputfile.open("../quad_info/quad_info_1d_p38.dat");
            break;
        case 1: inputfile.open("../quad_info/quad_info_2d_p20.dat");
            break;
        case 2: inputfile.open("../quad_info/nme6313-sup-0001-supinfo/NME_6313_cubature_tetra_p20_n448.dat");
            break;
        }
        for (int p = 0; p < n_q_point[ind]; ++p){
            for (int indt = 0; indt <= ind; ++indt)
                if (ind == DIM-1)
                    inputfile >> QPoint[p][indt];
                else
                    inputfile >> QPoint_Barycenter[ind][p][indt];
            if (ind == 0) // modify coordinate to be barycenter one
                QPoint_Barycenter[ind][p][0] = (QPoint_Barycenter[ind][p][0] + 1) * 0.5;
            if (ind < DIM-1){ // calculate the barycenter coordinate of the last point
                QPoint_Barycenter[ind][p][ind+1] = 1.0;
                for (int indt = 0; indt <= ind; ++indt)
                    QPoint_Barycenter[ind][p][ind+1] -= QPoint_Barycenter[ind][p][indt];
            }
            inputfile >> Weight[ind][p];
        }
        std::cout << "Read " << ind+1 << "d quadrature info with accuracy " << n_quad_accuracy[ind] << ", find " << n_q_point[ind] << " points and weights\n";
        inputfile.close();
    }

    
    // calculate the value of basis function at each quadrature point
    std::vector<std::vector<std::vector<double> > > basis_value; // [ind = 0:2][p = 0:n_q_point[ind]-1][]: basis_value for ind+1 dimensional quadrature region
    basis_value.resize(DIM);
    for (int ind = 0; ind < DIM; ++ind)
        basis_value[ind].resize(n_q_point[ind]);
    // basis function value at 3-d quadrature points
    for (int p = 0; p < n_q_point[DIM-1]; ++p){
        double x = QPoint[p][0], y = QPoint[p][1], z = QPoint[p][2];
        double xi = 2*x/(1-y-z)-1, eta = 2*y/(1-z)-1, zeta = 2*z-1;
        basis_value[DIM-1][p].resize(n_index);
        // calculate $J_k^{-1,-1}(xi)$ for k = 0: M
        double Jxi[M+1], Jeta[M+1], Jzeta[M+1];
        Jxi[0] = Jeta[0] = Jzeta[0] = 1;
        Jxi[1] = xi; // as J_1^{-1,-1}(xi) = xi
        for (int k = 1; k < M; ++k)
            Jxi[k+1] = ((xi - calc_coefficient_a(-1,-1,2,k)) * Jxi[k] - calc_coefficient_a(-1,-1,3,k) * Jxi[k-1]) / calc_coefficient_a(-1,-1,1,k);
        // traverse first component of multiindex
        for (int l1 = 0; l1 <= M; ++l1){
            int aph2 = 2 * l1 - 1; // alpha for the generalized jacobi polynomial of eta
            // calculate value of generalized jacobi polynomial Jeta
            Jeta[1] = calc_generalized_jacobi_polynomial(aph2, -1, 1, eta);
            for (int k = 1; k < M-l1; ++k)
                Jeta[k+1] = ((eta - calc_coefficient_a(aph2,-1,2,k)) * Jeta[k] - calc_coefficient_a(aph2,-1,3,k) * Jeta[k-1]) / calc_coefficient_a(aph2,-1,1,k);
            // traverse second component
            for (int l2 = 0; l2 <= M-l1; ++l2){
                int aph3 = 2 * l1 + 2 * l2 - 1;
                Jzeta[1] = calc_generalized_jacobi_polynomial(aph3, -1, 1, zeta);
                for (int k = 1; k < M-l1-l2; ++k)
                    Jzeta[k+1] = ((zeta - calc_coefficient_a(aph3,-1,2,k)) * Jzeta[k] - calc_coefficient_a(aph3,-1,3,k) * Jzeta[k-1]) / calc_coefficient_a(aph3,-1,1,k);
                // traverse third component
                for (int l3 = 0; l3 <= M-l1-l2; ++l3){
                    Multiindex index_now = Unitary_Multiindex[0] * l1 + Unitary_Multiindex[1] * l2 + Unitary_Multiindex[2] * l3;
                    basis_value[DIM-1][p][correspondence.index2number(index_now)] = pow(1-y-z,l1)*Jxi[l1] * pow(1-z,l2)*Jeta[l2] * Jzeta[l3];
                }
            }
        }
    }
    for (int p = 0; p < n_q_point[1]; ++p){
        basis_value[1][p].resize(n_dof_geometry[2]);
        double x = QPoint_Barycenter[1][p][0], y = QPoint_Barycenter[1][p][1];
        double xi = 2*x/(1-y)-1, eta = 2*y-1;
        double Jxi[M+1], Jeta[M+1];
        Jxi[0] = Jeta[0] = 1;
        Jxi[1] = xi;
        for (int l = 1; l < M; ++l)
            Jxi[l+1] = ((xi - calc_coefficient_a(-1,-1,2,l))*Jxi[l] - calc_coefficient_a(-1,-1,3,l)*Jxi[l-1]) / calc_coefficient_a(-1,-1,1,l);
        for (int l1 = 2; l1 <= M; ++l1){
            int aph = 2 * l1 - 1;
            Jeta[1] = calc_generalized_jacobi_polynomial(aph, -1, 1, eta);
            for (int l = 1; l < M; ++l)
                Jeta[l+1] = ((eta - calc_coefficient_a(aph,-1,2,l))*Jeta[l] - calc_coefficient_a(aph,-1,3,l)*Jeta[l-1]) / calc_coefficient_a(aph,-1,1,l);
            for (int l2 = 1; l2 <= M-l1; ++l2){
                int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2 - 1;
                basis_value[1][p][ind_index] = 2 * pow(1-y, l1) * Jxi[l1] * Jeta[l2];
            }
        }
    }
    for (int p = 0; p < n_q_point[0]; ++p){
        basis_value[0][p].resize(n_dof_geometry[1]);
        double xi = QPoint_Barycenter[0][p][0] * 2 - 1;
        double Jxi[M+1];
        Jxi[0] = 1;
        Jxi[1] = xi;
        for (int l = 1; l < M; ++l)
            Jxi[l+1] = ((xi - calc_coefficient_a(-1,-1,2,l))*Jxi[l] - calc_coefficient_a(-1,-1,3,l)*Jxi[l-1]) / calc_coefficient_a(-1,-1,1,l);
        for (int l = 2; l <= M; ++l)
            basis_value[0][p][l-2] = 2 * Jxi[l];
    }
    std::cerr << "calculate function value of generalized jacobi polynomial at quadrature points\n";
    // calculate the value of generalized jacobi polynomial for the interpolation of function
    std::vector<std::vector<std::vector<double > > > basis_value_interp(DIM);
    for (int ind = 0; ind < DIM; ++ind){
        basis_value_interp[ind].resize(n_q_point[ind]);
        for (int p = 0; p < n_q_point[ind]; ++p)
            basis_value_interp[ind][p].resize(n_dof_geometry[ind+1]);
    }
    // for the interpolation of the dof on edge
    for (int p = 0; p < n_q_point[0]; ++p){
        double xp = 2 * QPoint_Barycenter[0][p][0] - 1;
        basis_value_interp[0][p][0] = 1;
        basis_value_interp[0][p][1] = calc_generalized_jacobi_polynomial(1, 1, 1, xp);
        for (int l = 1; l < M-2; ++l)
            basis_value_interp[0][p][l+1] = ((xp - calc_coefficient_a(1,1,2,l)) * basis_value_interp[0][p][l]
                                             - calc_coefficient_a(1,1,3,l) * basis_value_interp[0][p][l-1]) / calc_coefficient_a(1,1,1,l);
    }
    // for the interpolation of the dof on face
    // additional basis value, [p][l][0]: $(1-y)^l*J_{l-2}^{1,1}(\xi_p)$, [p][l][1]: $J_{l-2}^{1,1}(\eta_p)$
    std::vector<std::vector<std::vector<double> > > basis_value_addition(n_q_point[1]); 
    for (int p = 0; p < n_q_point[1]; ++p){
        basis_value_addition[p].resize(n_dof_geometry[1]);
        for (int l = 0; l < n_dof_geometry[1]; ++l)
            basis_value_addition[p][l].resize(2);
    }
    for (int p = 0; p < n_q_point[1]; ++p){
        double xp = QPoint_Barycenter[1][p][0], yp = QPoint_Barycenter[1][p][1];
        double xi = 2*xp/(1-yp) - 1, eta = 2*yp - 1;
        double Jxi[M-1], Jeta[M];
        // evaluate basis_value_interp
        Jxi[0] = Jeta[0] = 1;
        Jxi[1] = calc_generalized_jacobi_polynomial(1, 1, 1, xi);
        for (int l = 1; l < M-2; ++l) // evaluate Jxi by recursion relation
            Jxi[l+1] = ((xi - calc_coefficient_a(1,1,2,l)) * Jxi[l] - calc_coefficient_a(1,1,3,l) * Jxi[l-1]) / calc_coefficient_a(1,1,1,l);
        for (int l1 = 2; l1 <= M; ++l1){
            Jeta[1] = calc_generalized_jacobi_polynomial(2*l1-1, 1, 1, eta);
            for (int l2 = 1; l2 < M-l1; ++l2) // evaluate Jeta by recursion relation, Jeta[M-l1] is un-used
                Jeta[l2+1] = ((eta - calc_coefficient_a(2*l1-1,1,2,l2)) * Jeta[l2] - calc_coefficient_a(2*l1-1,1,3,l2) * Jeta[l2-1]) / calc_coefficient_a(2*l1-1,1,1,l2);
            for (int l2 = 1; l2 <= M-l1; ++l2){ // assign basis_value_interp, corresponds to (l1-2, l2-1)
                int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2-1;
                basis_value_interp[1][p][ind_index] = pow(1-eta, l1-2) * Jxi[l1-2] * Jeta[l2-1];
            }
        }
        // evaluate basis_value_addition, which correspond to $(1-yp)^{l-2}J_{l-2}^{1,1}(xi)$ and $J_{l-2}^{1,1}(eta)$ for l = 2: M
        for (int l = 2; l <= M; ++l)
            basis_value_addition[p][l-2][0] = pow(1-yp, l-2) * Jxi[l-2]; // as Jxi is the same
        Jeta[1] = calc_generalized_jacobi_polynomial(1, 1, 1, eta);
        for (int l = 1; l < M-2; ++l)
            Jeta[l+1] = ((eta - calc_coefficient_a(1,1,2,l)) * Jeta[l] - calc_coefficient_a(1,1,3,l) * Jeta[l-1]) / calc_coefficient_a(1,1,1,l);
        for (int l = 2; l <= M; ++l)
            basis_value_addition[p][l-2][1] = Jeta[l-2];
    }
    // for the interpolation of the interior dof
    for (int p = 0; p < n_q_point[2]; ++p){
        double xp = QPoint[p][0], yp = QPoint[p][1], zp = QPoint[p][2];
        double xi = 2*xp/(1-yp-zp)-1, eta = 2*yp/(1-zp)-1, zeta = 2*zp-1;
        double Jxi[M-1], Jeta[M], Jzeta[M];
        Jxi[0] = Jeta[0] = Jzeta[0] = 1;
        Jxi[1] = calc_generalized_jacobi_polynomial(1, 1, 1, xi);
        for (int l = 1; l < M-2; ++l) // l1 = 2: M -> l1-2 = 0: M-2
            Jxi[l+1] = ((xi - calc_coefficient_a(1,1,2,l)) * Jxi[l] - calc_coefficient_a(1,1,3,l) * Jxi[l-1]) / calc_coefficient_a(1,1,1,l);
        for (int l1 = 2; l1 <= M; ++l1){
            Jeta[1] = calc_generalized_jacobi_polynomial(2*l1-1, 1, 1, eta);
            for (int l = 1; l < M - l1; ++l) // in fact, we consider l2 = 1: M-l1 -> l2-1 = 0: M-l1-1, so Jeta[M-l1] is un-used
                Jeta[l+1] = ((eta - calc_coefficient_a(2*l1-1,1,2,l)) * Jeta[l] - calc_coefficient_a(2*l1-1,1,3,l) * Jeta[l-1]) / calc_coefficient_a(2*l1-1,1,1,l);
            for (int l2 = 1; l2 <= M - l1; ++l2){
                int aph = 2 * l1 + 2 * l2 - 1;
                Jzeta[1] = calc_generalized_jacobi_polynomial(aph, 1, 1, zeta);
                for (int l = 1; l < M - l1 - l2; ++l) // similarly, Jzeta[M-l1-l2] is un-used
                    Jzeta[l+1] = ((zeta - calc_coefficient_a(aph,1,2,l)) * Jzeta[l] - calc_coefficient_a(aph,1,3,l) * Jzeta[l-1]) / calc_coefficient_a(aph,1,1,l);
                for (int l3 = 1; l3 <= M - l1 - l2; ++l3){
                    int ind_index = correspondence.index2number(Unitary_Multiindex[0] * (l1-2) + Unitary_Multiindex[1] * (l2-1) + Unitary_Multiindex[2] * (l3-1));
                    basis_value_interp[2][p][ind_index] = Jxi[l1-2] * Jeta[l2-1]*pow(1-eta,l1-2) * Jzeta[l3-1]*pow(1-zeta,l1+l2-3);
                }
            }
        }
    }
    std::cerr << "calculate the basis value for interpolation\n";
    
    
    // assign the weight_location, geometry_dimension and geometry_order of dof, which determine the order of these geometry in discretized matrix
    int n_geometry_total = 0;
    for (int ind = 0; ind <= DIM; ++ind)
        n_geometry_total += n_geometry[ind];
    std::vector<double> weight_location(n_geometry_total); // the weight_location of 0: 3 dimensional geometry in turns
    std::vector<int> geometry_dimension(n_geometry_total);
    std::vector<int> geometry_order(n_geometry_total); // [i = 0:n_geometry_total-1] the order of i-th entry in weight_location, whose dimension is geometry_dimension[i]
    // correspondence between fem dof/element and geometry
    std::vector<int> transform_femdof2geometry(n_dof, -1); // index from fem dof to 0 & 1 dimensional geometry, for 1 dimensional geometry, its index plus mesh.n_geometry(0)
    std::vector<std::vector<double> > transform_femele2geometry(n_element); // [i][ind = 0:4] index of geometry in i-th fem element, ind = 0:3 correspond to face (in front of ind-th vertex), ind = 4 for tetrahedron
    for (int i = 0; i < n_element; ++i)
        transform_femele2geometry[i].resize(5);
    std::vector<int> location_geometry(n_geometry_total); // location of all geometry (0:3 dimensional) according to increasing order of weight_location
    std::vector<int> location_actualdof(n_geometry_total); // start index of geometry in actual discretized matrix
    /*
     * transform_femdof2geometry (0 & 1 dimensional geometry) -> total index = 0: n_geometry_total -> location_actualdof[location_geometry[total index]]
     * transform_femele2geometry (2 & 3 dimensional geometry)
     */
    for (int ind = 0; ind <= DIM; ++ind){
        int index_start = 0;
        for (int indt = 0; indt < ind; ++indt)
            index_start += n_geometry[indt];
        for (int i = 0; i < n_geometry[ind]; ++i){
            geometry_dimension[index_start + i] = ind;
            geometry_order[index_start + i] = i;
        }
    }
    // for (int i = 0; i < n_geometry_total; ++i)
    //     std::cout << "the " << i << "-th geometry in total order has dimension " << geometry_dimension[i] << ", and order " << geometry_order[i] <<  '\n';
    // construct correspondence between fem dof and geometry (I): find weight_location for 0 & 1 dimensional geometry according to the order in fem_space
    for (int i = 0; i < n_dof; ++i)
        for (int ind = 0; ind <= 1 && transform_femdof2geometry[i] < 0; ++ind) // as all order in fem_space is natural number, so >= 0
            for (int j = 0; j < n_geometry[ind]; ++j)
                if (distance(point_ref_mesh[ind][j], fem_space.dofInfo(i).interp_point) < 1.0e-8){
                    transform_femdof2geometry[i] = ind * n_geometry[0] + j; // if ind == 1, add n_geometry[0], total index of point or edge
                    weight_location[transform_femdof2geometry[i]] = i;
                    break;
                }
    // calculate weight_location for 2 & 3 dimensional geometry
    for (int ind = 2; ind <= DIM; ++ind){
        int index_start = n_geometry[0] + n_geometry[1] + (ind-2) * n_geometry[2]; // if ind == 3, add n_geometry[2]
        for (int i = 0; i < n_geometry[ind]; ++i){
            // std::cout << "the " << i << "-th " << ind << " dimensional geomtry\n";
            double weight_tmp = 0.0;
            for (int indt = 0; indt <= ind; ++indt){
                // std::cout << "\tthe " << indt << "-th vertex of this geomtry has index " << mesh.geometry(ind, i).vertex(indt)
                //           << " with weight_location " << weight_location[mesh.geometry(ind, i).vertex(indt)] << '\n';
                weight_tmp += weight_location[mesh.geometry(ind, i).vertex(indt)];
            }
            weight_tmp /= (ind + 1.0);
            weight_location[index_start + i] = weight_tmp;
            // std::cout << "\t\tnow weight_tmp = " << weight_tmp << '\n';
        }
    }
    // for (int i = 0; i < n_geometry_total; ++i)
    //     std::cout << "the " << i << "-th geometry has weight_location " << weight_location[i] << '\n';
    // sort according to weight_location
    for (int i = 1; i < n_geometry_total; ++i)
        for (int j = i; j > 0; --j)
            if (weight_location[j] < weight_location[j-1]){
                double tmp = weight_location[j];
                weight_location[j] = weight_location[j-1];
                weight_location[j-1] = tmp;
                int tmpi = geometry_dimension[j];
                geometry_dimension[j] = geometry_dimension[j-1];
                geometry_dimension[j-1] = tmpi;
                tmpi = geometry_order[j];
                geometry_order[j] = geometry_order[j-1];
                geometry_order[j-1] = tmpi;
            }
    for (int i = 1; i < n_geometry_total; ++i)
        if (weight_location[i] < weight_location[i-1])
            std::cout << "error in sort of weight_location\n";
    for (int i = 0; i < n_geometry_total; ++i){ // assign location for geometry_order[i]-th geometry_dimension[i] dimensional geometry
        int index = geometry_order[i]; // recover the index of geometry in whole order
        for (int ind = 0; ind < geometry_dimension[i]; ++ind)
            index += n_geometry[ind];
        location_geometry[index] = i;
    }
    // insert dof into each location_geometry
    location_actualdof[0] = 0;
    for (int i = 1; i < n_geometry_total; ++i)
        location_actualdof[i] = location_actualdof[i-1] + n_dof_geometry[geometry_dimension[i-1]];
    // construct correspondence between fem element and geometry (II), traverse fem element, find location for 2 & 3 dimensional geometry of fem element
    the_element = fem_space.beginElement();
    for (int i = 0; the_element != end_element; ++the_element, ++i){
        const std::vector<int>& element_dof = the_element->dof();
        AFEPack::Point<DIM> point_tmp;
        for (int j = 0; j < 4; ++j) // the summation of all points
            point_tmp += fem_space.dofInfo(element_dof[j]).interp_point;
        AFEPack::Point<DIM> point_now = point_tmp;
        point_now *= 0.25;
        int start_tetrahedron = n_geometry[0] + n_geometry[1] + n_geometry[2];
        for (int j = 0; j < n_geometry[3]; ++j)
            if (distance(point_now, point_ref_mesh[3][j]) < 1.0e-8){
                transform_femele2geometry[i][4] = start_tetrahedron + j; // total index of tetrahedron
                break;
            }
        int start_triangle = n_geometry[0] + n_geometry[1];
        for (int ind = 0; ind < 4; ++ind){
            point_now = point_tmp;
            point_now -= fem_space.dofInfo(element_dof[ind]).interp_point;
            point_now /= 3.0;
            for (int j = 0; j < n_geometry[2]; ++j)
                if (distance(point_now, point_ref_mesh[2][j]) < 1.0e-8){
                    transform_femele2geometry[i][ind] = start_triangle + j; // total index of triangle
                    break;
                }
        }
    }

    // construct index for dof in each fem element
    std::vector<std::vector<double> > transform_fem2dof(n_element);
    for (int i = 0; i < n_element; ++i)
        transform_fem2dof[i].resize(n_index);
    the_element = fem_space.beginElement();
    for (int i = 0; the_element != end_element; ++the_element, ++i){
        const std::vector<int>& element_dof = the_element->dof();
        // point index, (0, 0, 0) and Unitary_Multiindex[0:2]
        for (int j = 0; j <= 3; ++j)
            transform_fem2dof[i][j] = location_actualdof[location_geometry[transform_femdof2geometry[element_dof[j]]]];
        // edge index
        for (int ind = 0; ind < DIM; ++ind){
            // $E_{01}$, $E_{02}$, $E_{03}$, correspond to element_dof[4:6], respectively
            int location_start = location_actualdof[location_geometry[transform_femdof2geometry[element_dof[ind+4]]]];
            for (int l = 2; l <= M; ++l)
                transform_fem2dof[i][correspondence.index2number(Unitary_Multiindex[ind] * l)] = location_start + l - 2;
            // $E_{23}$, $E_{13}$, $E_{12}$, correspond to element_dof[7:9], respectively
            location_start = location_actualdof[location_geometry[transform_femdof2geometry[element_dof[ind+7]]]];
            int index_smaller = ((ind == 0) ? 2 : 1); // smaller nonzero index of multiindex: 23 -> 2, 13 -> 1, 12 -> 1
            for (int l = 2; l <= M; ++l){
                Multiindex index_tmp = Unitary_Multiindex[index_smaller -1] + Unitary_Multiindex[5-ind - index_smaller -1] * (l-1); // 5-ind = summation of indices of two endpoints, -1 for Unitary_Multiindex[0:2] but not 1:3
                transform_fem2dof[i][correspondence.index2number(index_tmp)] = location_start + l - 2;
            }
        }
        // face index
        for (int ind = 0; ind <= 3; ++ind){
            int location_start = location_actualdof[location_geometry[transform_femele2geometry[i][ind]]];
            for (int l1 = 2; l1 <= M; ++l1)
                for (int l2 = 1; l2 <= M-l1; ++l2){
                    int order = (1+l1+l2-3)*(l1+l2-3)/2 + l2-1; // order of multiindex (l1, l2) with l1 >= 2, l2 >= 1, which is equivalent to the lexicographical order of (l1-2, l2-1)
                    Multiindex index_tmp;
                    if (ind == 0)
                        index_tmp = Unitary_Multiindex[0] + Unitary_Multiindex[1] * (l1-1) + Unitary_Multiindex[2] * l2;
                    else{
                        int index_smaller = ((ind == 1) ? 2 : 1); // smaller location nonzero index of multiindex
                        int index_larger  = ((ind == 3) ? 2 : 3); // larger location nonzero index of multiindex
                        index_tmp = Unitary_Multiindex[index_smaller-1] * l1 + Unitary_Multiindex[index_larger-1] * l2;
                    }
                    transform_fem2dof[i][correspondence.index2number(index_tmp)] = location_start + order;
                }
        }
        // tetrahedron index
        int location_start_tetrahedron = location_actualdof[location_geometry[transform_femele2geometry[i][4]]];
        for (int l1 = 2; l1 <= M; ++l1)
            for (int l2 = 1; l2 <= M; ++l2)
                for (int l3 = 1; l3 <= M; ++l3){
                    if (l1 + l2 + l3 > M) continue;
                    Multiindex index_tmp = Unitary_Multiindex[0] * l1 + Unitary_Multiindex[1] * l2 + Unitary_Multiindex[2] * l3;
                    int order = correspondence.index2number(index_tmp - Unitary_Multiindex[0] * 2 - Unitary_Multiindex[1] - Unitary_Multiindex[2]);
                    transform_fem2dof[i][correspondence.index2number(index_tmp)] = location_start_tetrahedron + order;
                }
    }

    // construct tranform_local, weight_transform_local: generalize jacobi polynomial -> basis function
    n_transform_local.resize(n_index);
    transform_local.resize(n_index);
    weight_transform_local.resize(n_index);
    // vertex modes
    n_transform_local[0] = 4; n_transform_local[1] = 2; n_transform_local[2] = 3; n_transform_local[3] = 4;
    for (int ind = 0; ind <= 3; ++ind){
        transform_local[ind].resize(n_transform_local[ind]);
        weight_transform_local[ind].resize(n_transform_local[ind]);
    }
    transform_local[0][0] = 0; weight_transform_local[0][0] = 0.125; // $J_{0,0,0} -> \varphi_{0,0,0}$
    transform_local[0][1] = 1; weight_transform_local[0][1] = 0.125; // $J_{0,0,0} -> \varphi_{1,0,0}$
    transform_local[0][2] = 2; weight_transform_local[0][2] = 0.25;  // $J_{0,0,0} -> \varphi_{0,1,0}$
    transform_local[0][3] = 3; weight_transform_local[0][3] = 0.5;   // $J_{0,0,0} -> \varphi_{0,0,1}$
    transform_local[1][0] = 0; weight_transform_local[1][0] = -0.5;  // $J_{1,0,0} -> \varphi_{0,0,0}$
    transform_local[1][1] = 1; weight_transform_local[1][1] = 0.5;   // $J_{1,0,0} -> \varphi_{1,0,0}$
    transform_local[2][0] = 0; weight_transform_local[2][0] = -0.25; // $J_{0,1,0} -> \varphi_{0,0,0}$
    transform_local[2][1] = 1; weight_transform_local[2][1] = -0.25; // $J_{0,1,0} -> \varphi_{1,0,0}$
    transform_local[2][2] = 2; weight_transform_local[2][2] = 0.5;   // $J_{1,1,0} -> \varphi_{0,1,0}$
    transform_local[3][0] = 0; weight_transform_local[3][0] = -0.125;// $J_{0,0,1} -> \varphi_{0,0,0}$
    transform_local[3][1] = 1; weight_transform_local[3][1] = -0.125;// $J_{0,0,1} -> \varphi_{1,0,0}$
    transform_local[3][2] = 2; weight_transform_local[3][2] = -0.25; // $J_{0,0,1} -> \varphi_{0,1,0}$
    transform_local[3][3] = 3; weight_transform_local[3][3] = 0.5;   // $J_{0,0,1} -> \varphi_{0,0,1}$
    // edge modes
    for (int l = 2; l <= M; ++l){
        // $J_{l,0,0}$ -> $\varphi_{l,0,0}$ 
        Multiindex index_tmp3 = Unitary_Multiindex[0] * l;
        int ind_tmp3 = correspondence.index2number(index_tmp3);
        n_transform_local[ind_tmp3] = 1;
        transform_local[ind_tmp3].resize(n_transform_local[ind_tmp3]);
        weight_transform_local[ind_tmp3].resize(n_transform_local[ind_tmp3]);
        transform_local[ind_tmp3][0] = ind_tmp3; weight_transform_local[ind_tmp3][0] = 2;
        // $J_{0,l,0}$, $J_{1,l-1,0}$ -> $\varphi_{0,l,0}$, $\varphi_{1,l-1,0}$
        Multiindex index_tmp1 = Unitary_Multiindex[1] * l;
        Multiindex index_tmp2 = Unitary_Multiindex[0] + Unitary_Multiindex[1] * (l-1);
        int ind_tmp1 = correspondence.index2number(index_tmp1);
        int ind_tmp2 = correspondence.index2number(index_tmp2);
        n_transform_local[ind_tmp1] = 2;
        n_transform_local[ind_tmp2] = 2;
        transform_local[ind_tmp1].resize(n_transform_local[ind_tmp1]);
        transform_local[ind_tmp2].resize(n_transform_local[ind_tmp2]);
        weight_transform_local[ind_tmp1].resize(n_transform_local[ind_tmp1]);
        weight_transform_local[ind_tmp2].resize(n_transform_local[ind_tmp2]);
        transform_local[ind_tmp1][0] = ind_tmp1; weight_transform_local[ind_tmp1][0] = 1;         // $J_{0,l,0}$   -> $\varphi_{0,l,0}$
        transform_local[ind_tmp1][1] = ind_tmp2; weight_transform_local[ind_tmp1][1] = 1;         // $J_{0,l,0}$   -> $\varphi_{1,l-1,0}$
        transform_local[ind_tmp2][0] = ind_tmp1; weight_transform_local[ind_tmp2][0] = (l-1.0)/l; // $J_{1,l-1,0}$ -> $\varphi_{0,l,0}$
        transform_local[ind_tmp2][1] = ind_tmp2; weight_transform_local[ind_tmp2][1] =-(l-1.0)/l; // $J_{1,l-1,0}$ -> $\varphi_{1,l-1,0}$
        // $J_{0,0,l}$, $J_{1,0,l-1}$, $J_{0,1,l-1}$ -> $\varphi_{0,0,l}$, $varphi_{1,0,l-1}$, $varphi_{0,1,l-1}$
        index_tmp1 = Unitary_Multiindex[2] * l;
        index_tmp2 = Unitary_Multiindex[0] + Unitary_Multiindex[2] * (l-1);
        index_tmp3 = Unitary_Multiindex[1] + Unitary_Multiindex[2] * (l-1);
        ind_tmp1 = correspondence.index2number(index_tmp1);
        ind_tmp2 = correspondence.index2number(index_tmp2);
        ind_tmp3 = correspondence.index2number(index_tmp3);
        n_transform_local[ind_tmp1] = 3;
        n_transform_local[ind_tmp2] = 2;
        n_transform_local[ind_tmp3] = 3;
        transform_local[ind_tmp1].resize(n_transform_local[ind_tmp1]);
        transform_local[ind_tmp2].resize(n_transform_local[ind_tmp2]);
        transform_local[ind_tmp3].resize(n_transform_local[ind_tmp3]);
        weight_transform_local[ind_tmp1].resize(n_transform_local[ind_tmp1]);
        weight_transform_local[ind_tmp2].resize(n_transform_local[ind_tmp2]);
        weight_transform_local[ind_tmp3].resize(n_transform_local[ind_tmp3]);
        transform_local[ind_tmp1][0] = ind_tmp1; weight_transform_local[ind_tmp1][0] = 0.5;           // $J_{0,0,l}$   -> $\varphi_{0,0,l}$
        transform_local[ind_tmp1][1] = ind_tmp2; weight_transform_local[ind_tmp1][1] = 0.5;           // $J_{0,0,l}$   -> $\varphi_{1,0,l-1}$
        transform_local[ind_tmp1][2] = ind_tmp3; weight_transform_local[ind_tmp1][2] = 1;             // $J_{0,0,l}$   -> $\varphi_{0,1,l-1}$
        transform_local[ind_tmp2][0] = ind_tmp1; weight_transform_local[ind_tmp2][0] = (l-1.0)/l;     // $J_{1,0,l-1}$ -> $\varphi_{0,0,l}$
        transform_local[ind_tmp2][1] = ind_tmp2; weight_transform_local[ind_tmp2][1] =-(l-1.0)/l;     // $J_{1,0,l-1}$ -> $\varphi_{1,0,l-1}$
        transform_local[ind_tmp3][0] = ind_tmp1; weight_transform_local[ind_tmp3][0] = (l-1.0)/(2*l); // $J_{0,1,l-1}$ -> $\varphi_{0,0,l}$
        transform_local[ind_tmp3][1] = ind_tmp2; weight_transform_local[ind_tmp3][1] = (l-1.0)/(2*l); // $J_{0,1,l-1}$ -> $\varphi_{1,0,l-1}$
        transform_local[ind_tmp3][2] = ind_tmp3; weight_transform_local[ind_tmp3][2] =-(l-1.0)/l;     // $J_{0,1,l-1}$ -> $\varphi_{0,1,l-1}$
    }
    // face modes
    for (int l1 = 2; l1 <= M; ++l1)
        for (int l2 = 1; l2 <= M-l1; ++l2){
            // $J_{0,l1,l2}$, $J_{1,l1-1,l2}$ -> $\varphi_{0,l1,l2}$, $varphi_{1,l1-1,l2}$
            Multiindex index_tmp1 = Unitary_Multiindex[1] * l1 + Unitary_Multiindex[2] * l2;
            Multiindex index_tmp2 = Unitary_Multiindex[0] + Unitary_Multiindex[1] * (l1-1) + Unitary_Multiindex[2] * l2;
            int ind_tmp1 = correspondence.index2number(index_tmp1);
            int ind_tmp2 = correspondence.index2number(index_tmp2);
            n_transform_local[ind_tmp1] = 2;
            n_transform_local[ind_tmp2] = 2;
            transform_local[ind_tmp1].resize(n_transform_local[ind_tmp1]);
            transform_local[ind_tmp2].resize(n_transform_local[ind_tmp2]);
            weight_transform_local[ind_tmp1].resize(n_transform_local[ind_tmp1]);
            weight_transform_local[ind_tmp2].resize(n_transform_local[ind_tmp2]);
            transform_local[ind_tmp1][0] = ind_tmp1; weight_transform_local[ind_tmp1][0] = 1;           // $J_{0,l1,l2}$   -> $\varphi_{0,l1,l2}$
            transform_local[ind_tmp1][1] = ind_tmp2; weight_transform_local[ind_tmp1][1] = 1;           // $J_{0,l1,l2}$   -> $\varphi_{1,l1-1,l2}$
            transform_local[ind_tmp2][0] = ind_tmp1; weight_transform_local[ind_tmp2][0] = (l1-1.0)/l1; // $J_{1,l1-1,l2}$ -> $\varphi_{0,l1,l2}$
            transform_local[ind_tmp2][1] = ind_tmp2; weight_transform_local[ind_tmp2][1] =-(l1-1.0)/l1; // $J_{1,l1-1,l2}$ -> $\varphi_{1,l1-1,l2}$
            // $J_{l1,0,l2}$ -> $\varphi_{l1,0,l2}$ and $J_{l1,l2,0}$ -> $\varphi_{l1,l2,0}$
            for (int ind = 1; ind <= 2; ++ind){
                Multiindex index_tmp = Unitary_Multiindex[0] * l1 + Unitary_Multiindex[ind] * l2;
                int ind_tmp = correspondence.index2number(index_tmp);
                n_transform_local[ind_tmp] = 1;
                transform_local[ind_tmp].resize(n_transform_local[ind_tmp]);
                weight_transform_local[ind_tmp].resize(n_transform_local[ind_tmp]);
                transform_local[ind_tmp][0] = ind_tmp; weight_transform_local[ind_tmp][0] = 2;
            }
        }
    // interior modes
    for (int l1 = 2; l1 <= M; ++l1)
        for (int l2 = 1; l2 <= M-l1; ++l2)
            for (int l3 = 1; l3 <= M-l1-l2; ++l3){
                Multiindex index_tmp = Unitary_Multiindex[0] * l1 + Unitary_Multiindex[1] * l2 + Unitary_Multiindex[2] * l3;
                int ind_tmp = correspondence.index2number(index_tmp);
                n_transform_local[ind_tmp] = 1;
                transform_local[ind_tmp].resize(n_transform_local[ind_tmp]);
                weight_transform_local[ind_tmp].resize(n_transform_local[ind_tmp]);
                transform_local[ind_tmp][0] = ind_tmp; weight_transform_local[ind_tmp][0] = 1;
            }

    // generate actual basis function value at local fem element, by transform_local and weigth_transform_local
    std::vector<std::vector<double> > basis_value_actual(n_q_point[DIM-1]); // [p = 0:n_q_point[ind]-1][]: basis_value for 3 dimensional quadrature region
    for (int p = 0; p < n_q_point[DIM-1]; ++p){
        basis_value_actual[p].resize(n_index);
        for (int i = 0; i < n_index; ++i)
            basis_value_actual[p][i] = 0;
    }
    // traverse 3-d quadrature point, use transform_local and weight_transform_local, calculate basis_value_actual
    for (int p = 0; p < n_q_point[DIM-1]; ++p)
        for (int i = 0; i < n_index; ++i)
            for (int j = 0; j < n_transform_local[i]; ++j)
                basis_value_actual[p][transform_local[i][j]] += weight_transform_local[i][j] * basis_value[DIM-1][p][i];


    // generate linear finite element solution by h-adaptive method, with sem solution given by argv[3]
    IrregularMesh<DIM> irregular_mesh_output(h_tree);
    if (N_GlobalRefine > 0)
        irregular_mesh_output.globalRefine(N_GlobalRefine); // global refine
    
    std::vector<double> err_indicator(20000000); // store error indicator, l2 error of each element
    double rate_refine = 0.5;
    double convergenceCoefficient = pow(2, DIM), refine_threshold = 1.33333;
    std::vector<double> val_interp(5000000);

    double tol_refine_l2 = atof(argv[4]);
    int n_local_refine = atoi(argv[3]); // maximum number of local refine
    double err_max_local_l2 = 0, err_global_l2 = 0;
    for (int n_refine = 0; n_refine <= n_local_refine; ++n_refine){
        // except inital step, local refinement according to error indicator
        if (n_refine > 0){
            // copy error indicator
            RegularMesh<DIM>& regular_mesh_t = irregular_mesh_output.regularMesh();
            Indicator<DIM> indicator(regular_mesh_t);
            for (int i = 0; i < regular_mesh_t.n_geometry(DIM); ++i)
                indicator[i] = err_indicator[i];
            std::cerr << "read indicator from err_indicator last step\n";
            // mesh adaption
            MeshAdaptor<DIM> mesh_adaptor(irregular_mesh_output);
            mesh_adaptor.convergenceOrder() = 0.;
            mesh_adaptor.refineStep() = 0;
            mesh_adaptor.setIndicator(indicator);
            mesh_adaptor.tolerence() = err_max_local_l2 * rate_refine / (convergenceCoefficient * refine_threshold);
            mesh_adaptor.is_refine_only() = true;
            mesh_adaptor.adapt();
            std::cerr << "locally refine the current mesh\n";
        }
        // calculate the error indicator for current mesh
        irregular_mesh_output.semiregularize();
        irregular_mesh_output.regularize(false);
        RegularMesh<DIM>& regular_mesh = irregular_mesh_output.regularMesh();
        // interpolate the density to current mesh with linear finite element basis function
        // val_interp.resize(regular_mesh.n_geometry(0));
#pragma omp parallel for
        for (int i = 0; i < regular_mesh.n_geometry(0); ++i){
            val_interp[i] = 0;
            AFEPack::Point<DIM> p = regular_mesh.point(i);
            std::vector<AFEPack::Point<DIM> > vertex(4);
            bool flag = false; // flag whether val_interp[i] is assigned a value
            for (int j = 0; j < mesh.n_geometry(0); ++j){ // match vertex
                if (distance(p, mesh.point(j)) > 1.0e-8) continue;
                int location = location_actualdof[location_geometry[j]];
                val_interp[i] = sol_sem(location);
                flag = true;
                break;
            }
            if (flag) continue;
            for (int j = 0; j < mesh.n_geometry(1); ++j){ // match edge
                int ind_s = mesh.geometry(1, j).vertex(0);
                int ind_e = mesh.geometry(1, j).vertex(1);
                int loc_s = location_actualdof[location_geometry[ind_s]];
                int loc_e = location_actualdof[location_geometry[ind_e]];
                int location = location_actualdof[location_geometry[j + n_geometry[0]]];
                AFEPack::Point<DIM> ps = mesh.point(ind_s);
                AFEPack::Point<DIM> pe = mesh.point(ind_e);
                double d  = distance(ps, pe), ds = distance(ps, p), de = distance(pe, p);
                if (fabs(d - (ds+de)) > 1.0e-8) continue;
                double x = ds / d, r = 1 - x;
                double xi = 2 * x - 1;
                double Jxi[M+1];
                Jxi[0] = 1; Jxi[1] = xi;
                for (int l1 = 1; l1 < M; ++l1)
                    Jxi[l1+1]  = ((xi - calc_coefficient_a(-1, -1, 2, l1)) *  Jxi[l1] - calc_coefficient_a(-1, -1, 3, l1) *  Jxi[l1-1])
                        / calc_coefficient_a(-1, -1, 1, l1);
                val_interp[i] = sol_sem(loc_s) * r + sol_sem(loc_e) * x;
                for (int l1 = 2; l1 <= M; ++l1)
                    val_interp[i] += Jxi[l1] * sol_sem(location + l1-2);
                flag = true;
                break;
            }
            if (flag) continue;
            for (int j = 0; j < mesh.n_geometry(2); ++j){ // match face
                for (int ind = 0; ind < 3; ++ind)
                    vertex[ind] = mesh.point(mesh.geometry(2, j).vertex(ind));
                double volume  = calc_volume_triangle(vertex[0], vertex[1], vertex[2]);
                double volume1 = calc_volume_triangle(vertex[0], p,         vertex[2]);
                double volume2 = calc_volume_triangle(vertex[0], vertex[1], p);
                double volume3 = calc_volume_triangle(p,         vertex[1], vertex[2]);
                if (fabs(volume - (volume1+volume2+volume3)) > 1.0e-8) continue;
                int ind_v[3], ind_e[3]; // location of first dof on corresponding geometry
                for (int ind = 0; ind < 3; ++ind){
                    ind_v[ind] = location_actualdof[location_geometry[mesh.geometry(2, j).vertex(ind)]];
                    ind_e[ind] = location_actualdof[location_geometry[mesh.geometry(2, j).boundary(ind) + n_geometry[0]]];
                }
                int location = location_actualdof[location_geometry[j + n_geometry[0] + n_geometry[1]]];
                std::vector<double> bas_val(n_dof_geometry[2], 0);
                double x = volume1 / volume, y = volume2 / volume, r = 1-x-y;
                double xi = 2*x/(1-y)-1, eta = 2*y-1;
                double Jxi[M+1], Jeta[M+1];
                Jxi[0] = Jeta[0] = 1; Jxi[1] = xi;
                for (int l1 = 1; l1 < M; ++l1)
                    Jxi[l1+1]  = ((xi - calc_coefficient_a(-1, -1, 2, l1)) *  Jxi[l1] - calc_coefficient_a(-1, -1, 3, l1) *  Jxi[l1-1])
                        / calc_coefficient_a(-1, -1, 1, l1);
                for (int l1 = 2; l1 <= M-1; ++l1){
                    int aph2 = 2 * l1 - 1;
                    Jeta[1]  = calc_generalized_jacobi_polynomial( aph2, -1, 1, eta);
                    for (int l2 = 1; l2 < M-l1; ++l2)
                        Jeta[l2+1]  = ((eta - calc_coefficient_a( aph2, -1, 2, l2)) *  Jeta[l2] - calc_coefficient_a( aph2, -1, 3, l2) *  Jeta[l2-1])
                            / calc_coefficient_a( aph2, -1, 1, l2);
                    for (int l2 = 1; l2 <= M-l1; ++l2){
                        int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2-1;
                        bas_val[ind_index] = 2 * pow(1-y, l1) * Jxi[l1] * Jeta[l2];
                    }
                }
                std::vector<double> bas_val_add_xi( n_dof_geometry[1], 0);
                std::vector<double> bas_val_add_eta(n_dof_geometry[1], 0);
                Jxi[0] = Jeta[0] = 1;
                Jxi[1]  = calc_generalized_jacobi_polynomial(1, 1, 1, xi);
                Jeta[1] = calc_generalized_jacobi_polynomial(1, 1, 1, eta);
                for (int l = 1; l < M; ++l){
                    Jxi[ l+1] = ((xi  - calc_coefficient_a(1, 1, 2, l)) *  Jxi[l] - calc_coefficient_a(1, 1, 3, l) *  Jxi[l-1]) / calc_coefficient_a(1, 1, 1, l);
                    Jeta[l+1] = ((eta - calc_coefficient_a(1, 1, 2, l)) * Jeta[l] - calc_coefficient_a(1, 1, 3, l) * Jeta[l-1]) / calc_coefficient_a(1, 1, 1, l);
                }
                for (int l = 0; l < n_dof_geometry[1]; ++l){
                    bas_val_add_xi[l]  = Jxi[l] * pow(1-y, l);
                    bas_val_add_eta[l] = Jeta[l];
                }
                
                val_interp[i] = sol_sem(ind_v[0]) * r + sol_sem(ind_v[1]) * x + sol_sem(ind_v[2]) * y;
                for (int l = 0; l < n_dof_geometry[1]; ++l)
                    val_interp[i] -= 2 * (x*y * sol_sem(ind_e[0]+l) * bas_val_add_eta[l] +
                                          y*r * sol_sem(ind_e[1]+l) * bas_val_add_eta[l] +
                                          x*r * sol_sem(ind_e[2]+l) * bas_val_add_xi[l]);
                for (int l1 = 2; l1 <= M; ++l1)
                    for (int l2 = 1; l2 <= M-l1; ++l2){
                        int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2-1;
                        val_interp[i] += bas_val[ind_index] * sol_sem(location + ind_index);
                    }
                flag = true;
                break;
            }
            if (flag) continue;
            double volume;
            for (int j = 0; j < n_element; ++j){
                for (int ind_vertex = 0; ind_vertex < 4; ++ind_vertex)
                    vertex[ind_vertex] = mesh.point(mesh.geometry(DIM, j).vertex(ind_vertex));
                volume = calc_volume_tetrahedron(vertex[0], vertex[1], vertex[2], vertex[3]);
                double volume0 = calc_volume_tetrahedron(p, vertex[1], vertex[2], vertex[3]);
                double volume1 = calc_volume_tetrahedron(vertex[0], p, vertex[2], vertex[3]);
                double volume2 = calc_volume_tetrahedron(vertex[0], vertex[1], p, vertex[3]);
                double volume3 = calc_volume_tetrahedron(vertex[0], vertex[1], vertex[2], p);
                double volume_local = volume0 + volume1 + volume2 + volume3;
                if (fabs(volume_local-volume) > 1.0e-8) continue;
                double x = volume1 / volume, y = volume2 / volume, z = volume3 / volume;
                double xi = 2*x/(1-y-z)-1, eta = 2*y/(1-z)-1, zeta = 2*z-1;
                std::vector<double> sol_(n_index, 0);
                // calculate immediate variable
                double Jxi[M+1], Jeta[M+1], Jzeta[M+1];
                Jxi[0] = Jeta[0] = Jzeta[0] = 1;
                Jxi[1]  = xi;
                for (int l1 = 1; l1 < M; ++l1)
                    Jxi[l1+1]  = ((xi - calc_coefficient_a(-1, -1, 2, l1)) *  Jxi[l1] - calc_coefficient_a(-1, -1, 3, l1) *  Jxi[l1-1])
                        / calc_coefficient_a(-1, -1, 1, l1);
                for (int l1 = 0; l1 <= M; ++l1){
                    int aph2 = 2 * l1 - 1;
                    Jeta[1]  = calc_generalized_jacobi_polynomial( aph2, -1, 1, eta);
                    for (int l2 = 1; l2 < M-l1; ++l2)
                        Jeta[l2+1]  = ((eta - calc_coefficient_a( aph2, -1, 2, l2)) *  Jeta[l2] - calc_coefficient_a( aph2, -1, 3, l2) *  Jeta[l2-1])
                            / calc_coefficient_a( aph2, -1, 1, l2);
                    for (int l2 = 0; l2 <= M-l1; ++l2){
                        int aph3 = 2 * l1 + 2 * l2 - 1;
                        Jzeta[1]  = calc_generalized_jacobi_polynomial( aph3, -1, 1, zeta);
                        for (int l3 = 1; l3 < M-l1-l2; ++l3)
                            Jzeta[l3+1]  = ((zeta - calc_coefficient_a( aph3, -1, 2, l3)) *  Jzeta[l3] - calc_coefficient_a( aph3, -1, 3, l3) *  Jzeta[l3-1])
                                / calc_coefficient_a( aph3, -1, 1, l3);
                        for (int l3 = 0; l3 <= M-l1-l2; ++l3){
                            Multiindex index_now = Unitary_Multiindex[0] * l1 + Unitary_Multiindex[1] * l2 + Unitary_Multiindex[2] * l3;
                            int ind_index = correspondence.index2number(index_now);
                            sol_[ind_index] = pow(1-y-z, l1) * Jxi[l1] * pow(1-z, l2) * Jeta[l2] * Jzeta[l3];
                        }
                    }
                }
                std::vector<double> sol_actual(n_index, 0);
                for (int ind_index = 0; ind_index < n_index; ++ind_index)
                    for (int ind_tl = 0; ind_tl < n_transform_local[ind_index]; ++ind_tl)
                        sol_actual[transform_local[ind_index][ind_tl]] += weight_transform_local[ind_index][ind_tl] * sol_[ind_index];

                val_interp[i] = 0;
                for (int ind_index = 0; ind_index < n_index; ++ind_index)
                    val_interp[i] += sol_sem(transform_fem2dof[j][ind_index]) * sol_actual[ind_index];
                break;
            }
        }
        std::cerr << "interpolate sem solution\n";

        if (n_refine == n_local_refine) break;
        
        // update variable
        // err_indicator.resize(regular_mesh.n_geometry(DIM), 0);
        for (int i = 0; i < regular_mesh.n_geometry(DIM); ++i)
            err_indicator[i] = 0;
        
        // calculate err_indicator
        int n_3d_geometry = regular_mesh.n_geometry(3);
#pragma omp parallel for
        for (int i = 0; i < n_3d_geometry; ++i){
            std::vector<AFEPack::Point<DIM> > vertex_p(4); // vertex for parent
            std::vector<double> volume_p(5); // volume used in calculation of parent index
            // search parent of this tetrahedron
            int parent = -1;
            AFEPack::Point<DIM> p_barycenter;
            p_barycenter[0] = p_barycenter[1] = p_barycenter[2] = 0;
            for (int j = 0; j <= 1; ++j)
                p_barycenter += regular_mesh.point(regular_mesh.geometry(3, i).vertex(j));
            if (regular_mesh.geometry(3, i).n_vertex() == 5)
                for (int j = 3; j <= 4; ++j)
                    p_barycenter += regular_mesh.point(regular_mesh.geometry(3, i).vertex(j));
            else
                for (int j = 2; j <= 3; ++j)
                    p_barycenter += regular_mesh.point(regular_mesh.geometry(3, i).vertex(j));
            p_barycenter *= 0.25;
            for (int j = 0; j < mesh.n_geometry(3); ++j){
                for (int k = 0; k < 4; ++k)
                    vertex_p[k] = mesh.point(mesh.geometry(3, j).vertex(k));
                volume_p[0] = calc_volume_tetrahedron(vertex_p[0],  vertex_p[1],  vertex_p[2],  vertex_p[3]); // volume of this tetrahedron
                volume_p[1] = calc_volume_tetrahedron(vertex_p[0],  p_barycenter, vertex_p[2],  vertex_p[3]);
                volume_p[2] = calc_volume_tetrahedron(vertex_p[0],  vertex_p[1],  p_barycenter, vertex_p[3]);
                volume_p[3] = calc_volume_tetrahedron(vertex_p[0],  vertex_p[1],  vertex_p[2],  p_barycenter);
                volume_p[4] = calc_volume_tetrahedron(p_barycenter, vertex_p[1],  vertex_p[2],  vertex_p[3]);
                double volume_local = 0;
                for (int k = 1; k <= 4; ++k) volume_local += volume_p[k];
                if (fabs(volume_p[0] - volume_local) > 1.0e-8) continue;
                parent = j;
                break;
            }
            
            // copy local info of vertex and solution
            std::vector<AFEPack::Point<DIM> > vertex_local(4);
            std::vector<double> sol_local(4);
            std::vector<double> sol_sem_local(n_index);
            for (int ind_v = 0; ind_v < 4; ++ind_v){
                vertex_local[ind_v] = regular_mesh.point(regular_mesh.geometry(3, i).vertex(ind_v));
                sol_local[ind_v] = val_interp[regular_mesh.geometry(3, i).vertex(ind_v)];
            }
            for (int ind_index = 0; ind_index < n_index; ++ind_index)
                sol_sem_local[ind_index] = sol_sem(transform_fem2dof[parent][ind_index]);
            
            // calculate error
            double err_local = 0;
            switch (regular_mesh.geometry(3, i).n_vertex()){
            case 4:
                err_local = calc_err(vertex_local, vertex_p, QPoint, Weight[2], sol_local, sol_sem_local);
                err_indicator[i] = sqrt(err_local);
                break;
            case 5: // twin tetrahedron
                vertex_local[3] = regular_mesh.point(regular_mesh.geometry(3, i).vertex(4));
                sol_local[3] = val_interp[regular_mesh.geometry(3, i).vertex(4)];
                err_local = calc_err(vertex_local, vertex_p, QPoint, Weight[2], sol_local, sol_sem_local);
                err_indicator[i] = err_local;
                vertex_local[1] = regular_mesh.point(regular_mesh.geometry(3, i).vertex(3));
                sol_local[1] = val_interp[regular_mesh.geometry(3, i).vertex(3)];
                err_local = calc_err(vertex_local, vertex_p, QPoint, Weight[2], sol_local, sol_sem_local);
                err_indicator[i] += err_local;
                err_indicator[i] = sqrt(err_indicator[i]);
                break;
            case 7: // four tetrahedron
                vertex_local[2] = regular_mesh.point(regular_mesh.geometry(3, i).vertex(5));
                vertex_local[3] = regular_mesh.point(regular_mesh.geometry(3, i).vertex(6));
                sol_local[2] = val_interp[regular_mesh.geometry(3, i).vertex(5)];
                sol_local[3] = val_interp[regular_mesh.geometry(3, i).vertex(6)];                
                err_local = calc_err(vertex_local, vertex_p, QPoint, Weight[2], sol_local, sol_sem_local);
                err_indicator[i] = err_local;
                vertex_local[1] = regular_mesh.point(regular_mesh.geometry(3, i).vertex(2));
                vertex_local[2] = regular_mesh.point(regular_mesh.geometry(3, i).vertex(4));
                sol_local[1] = val_interp[regular_mesh.geometry(3, i).vertex(2)];
                sol_local[2] = val_interp[regular_mesh.geometry(3, i).vertex(4)];
                err_local = calc_err(vertex_local, vertex_p, QPoint, Weight[2], sol_local, sol_sem_local);
                err_indicator[i] += err_local;
                vertex_local[1] = regular_mesh.point(regular_mesh.geometry(3, i).vertex(3));
                vertex_local[3] = regular_mesh.point(regular_mesh.geometry(3, i).vertex(5));
                sol_local[1] = val_interp[regular_mesh.geometry(3, i).vertex(3)];
                sol_local[3] = val_interp[regular_mesh.geometry(3, i).vertex(5)];
                err_local = calc_err(vertex_local, vertex_p, QPoint, Weight[2], sol_local, sol_sem_local);
                err_indicator[i] += err_local;
                vertex_local[1] = regular_mesh.point(regular_mesh.geometry(3, i).vertex(6));
                sol_local[1] = val_interp[regular_mesh.geometry(3, i).vertex(6)];
                err_local = calc_err(vertex_local, vertex_p, QPoint, Weight[2], sol_local, sol_sem_local);
                err_indicator[i] += err_local;
                err_indicator[i] = sqrt(err_indicator[i]);
            }
        }

        // evaluate err_global_l2
        err_global_l2 = 0;
        for (int i = 0; i < n_3d_geometry; ++i)
            err_global_l2 += pow(err_indicator[i], 2);
        err_global_l2 = sqrt(err_global_l2);
        
        // calculate error indicator
        err_max_local_l2 = 0;
        for (int i = 0; i < n_3d_geometry; ++i)
            if (err_indicator[i] > err_max_local_l2)
                err_max_local_l2 = err_indicator[i];
        
        // output error
        std::cout << "n_refine = " << n_refine << ", l2 error = " << err_global_l2 << '\n';
        std::cout << "\tn_dof = " << regular_mesh.n_geometry(0) << ", n_element = " << regular_mesh.n_geometry(3) << '\n';
        if (err_global_l2 < tol_refine_l2) break;
    }

    
    RegularMesh<DIM>& regular_mesh_output = irregular_mesh_output.regularMesh();
    // regular_mesh_output.writeOpenDXData("mesh_refine_output.dx");

    
    // generate fem_space_output
    TemplateGeometry<DIM> template_geometry_tetrahedron;
    CoordTransform<DIM, DIM> coord_transform_tetrahedron;
    TemplateDOF<DIM> template_dof_tetrahedron;
    BasisFunctionAdmin<double, DIM, DIM> basis_function_tetrahedron;
    template_geometry_tetrahedron.readData("tetrahedron.tmp_geo");
    coord_transform_tetrahedron.readData("tetrahedron.crd_trs");
    template_dof_tetrahedron.reinit(template_geometry_tetrahedron);
    template_dof_tetrahedron.readData("tetrahedron.1.tmp_dof");
    basis_function_tetrahedron.reinit(template_dof_tetrahedron);
    basis_function_tetrahedron.readData("tetrahedron.1.bas_fun");
    std::vector<TemplateElement<double, DIM, DIM> > template_element_output(3);
    template_element_output[0].reinit(template_geometry_tetrahedron,
                                      template_dof_tetrahedron,
                                      coord_transform_tetrahedron,
                                      basis_function_tetrahedron);
    TemplateGeometry<DIM> template_geometry_twin_tetrahedron;
    CoordTransform<DIM, DIM> coord_transform_twin_tetrahedron;
    TemplateDOF<DIM> template_dof_twin_tetrahedron;
    BasisFunctionAdmin<double, DIM, DIM> basis_function_twin_tetrahedron;
    template_geometry_twin_tetrahedron.readData("twin_tetrahedron.tmp_geo");
    coord_transform_twin_tetrahedron.readData("twin_tetrahedron.crd_trs");
    template_dof_twin_tetrahedron.reinit(template_geometry_twin_tetrahedron);
    template_dof_twin_tetrahedron.readData("twin_tetrahedron.1.tmp_dof");
    basis_function_twin_tetrahedron.reinit(template_dof_twin_tetrahedron);
    basis_function_twin_tetrahedron.readData("twin_tetrahedron.1.bas_fun");
    template_element_output[1].reinit(template_geometry_twin_tetrahedron,
                                      template_dof_twin_tetrahedron,
                                      coord_transform_twin_tetrahedron,
                                      basis_function_twin_tetrahedron);
    TemplateGeometry<DIM> template_geometry_four_tetrahedron;
    CoordTransform<DIM, DIM> coord_transform_four_tetrahedron;
    TemplateDOF<DIM> template_dof_four_tetrahedron;
    BasisFunctionAdmin<double, DIM, DIM> basis_function_four_tetrahedron;
    template_geometry_four_tetrahedron.readData("four_tetrahedron.tmp_geo");
    coord_transform_four_tetrahedron.readData("four_tetrahedron.crd_trs");
    template_dof_four_tetrahedron.reinit(template_geometry_four_tetrahedron);
    template_dof_four_tetrahedron.readData("four_tetrahedron.1.tmp_dof");
    basis_function_four_tetrahedron.reinit(template_dof_four_tetrahedron);
    basis_function_four_tetrahedron.readData("four_tetrahedron.1.bas_fun");
    template_element_output[2].reinit(template_geometry_four_tetrahedron,
                                      template_dof_four_tetrahedron,
                                      coord_transform_four_tetrahedron,
                                      basis_function_four_tetrahedron);

    // construct finite element space
    FEMSpace<double, DIM> fem_space_output(regular_mesh_output, template_element_output);
    int n_element_output = regular_mesh_output.n_geometry(DIM);
    fem_space_output.element().resize(n_element_output);
    for (int i = 0; i < n_element_output; ++i)
        switch (regular_mesh_output.geometry(DIM, i).n_vertex()){
        case 4:  fem_space_output.element(i).reinit(fem_space_output, i, 0);
            break;
        case 5:  fem_space_output.element(i).reinit(fem_space_output, i, 1);
            break;
        default: fem_space_output.element(i).reinit(fem_space_output, i, 2);
        } 
    fem_space_output.buildElement();
    fem_space_output.buildDof();
    fem_space_output.buildDofBoundaryMark();
    std::cout << "form fem_space for local refine mesh, find n_dof " << fem_space_output.n_dof() << '\n';

    FEMFunction<double, DIM> sol_output(fem_space_output);
    for (int i = 0; i < fem_space_output.n_dof(); ++i){
        AFEPack::Point<DIM> p_now = fem_space_output.dofInfo(i).interp_point;
        if (distance(p_now, regular_mesh_output.point(i)) < 1.0e-8){
            sol_output(i) = val_interp[i];
            continue;
        }
        for (int j = 0; j < regular_mesh_output.n_geometry(0); ++j)
            if (distance(p_now, regular_mesh_output.point(j)) < 1.0e-8){
                sol_output(i) = val_interp[j];
                break;
            }
    }
    sol_output.writeOpenDXData("sol_output.dx");
    
    return 0;
}

double calc_err(const std::vector<AFEPack::Point<DIM> >& vertex, const std::vector<AFEPack::Point<DIM> >& vertex_p,
                const std::vector<AFEPack::Point<DIM> >& QPoint, const std::vector<double>& Weight,
                const std::vector<double>& sol_local, const std::vector<double>& sol_sem_local)
{/* calculate l2 error
    evaluated with the aid of quadrature points on tetrahedron with nodes [vertex], node value [sol_local]
    parent tetrahedron with nodes [vertex_p], reference solution given by [sol_sem_local]
  */    
    int n_q_point = QPoint.size();
    int n_index = sol_sem_local.size();
    double err_local = 0;
    double err_local_l2 = 0;
    double volume_local = calc_volume_tetrahedron(vertex[0], vertex[1], vertex[2], vertex[3]);
    double volume = calc_volume_tetrahedron(vertex_p[0], vertex_p[1], vertex_p[2], vertex_p[3]);
    for (int i = 0; i < n_q_point; ++i){
        // find coordinates of this quadrature point, evaluate fem solution at this point
        AFEPack::Point<DIM> p;
        double sol_qpoint = 0;
        p[0] = p[1] = p[2] = 0;
        for (int ind_v = 0; ind_v < 4; ++ind_v){
            AFEPack::Point<DIM> p_tmp = vertex[ind_v];
            double rate = ind_v == 0 ? 1-QPoint[i][0]-QPoint[i][1]-QPoint[i][2] : QPoint[i][ind_v-1];
            p_tmp *= rate; p += p_tmp;
            sol_qpoint += sol_local[ind_v] * rate;
        }
        
        // evaluate sem solution at parent element
        double sol_sem_qpoint = 0;
        double x = calc_volume_tetrahedron(vertex_p[0], p, vertex_p[2], vertex_p[3]) / volume;
        double y = calc_volume_tetrahedron(vertex_p[0], vertex_p[1], p, vertex_p[3]) / volume;
        double z = calc_volume_tetrahedron(vertex_p[0], vertex_p[1], vertex_p[2], p) / volume;
        double xi = 2*x/(1-y-z)-1, eta = 2*y/(1-z)-1, zeta = 2*z-1;
        // calculate the contribution of inner product of gradient
        std::vector<double> sol_(n_index, 0);
        // calculate immediate variable
        double Jxi[M+1], Jeta[M+1], Jzeta[M+1];
        Jxi[0] = Jeta[0] = Jzeta[0] = 1;
        Jxi[1]  = xi;
        for (int l1 = 1; l1 < M; ++l1)
            Jxi[l1+1]  = ((xi - calc_coefficient_a(-1, -1, 2, l1)) *  Jxi[l1] - calc_coefficient_a(-1, -1, 3, l1) *  Jxi[l1-1])
                / calc_coefficient_a(-1, -1, 1, l1);
        for (int l1 = 0; l1 <= M; ++l1){
            int aph2 = 2 * l1 - 1;
            Jeta[1]  = calc_generalized_jacobi_polynomial( aph2, -1, 1, eta);
            for (int l2 = 1; l2 < M-l1; ++l2)
                Jeta[l2+1]  = ((eta - calc_coefficient_a( aph2, -1, 2, l2)) *  Jeta[l2] - calc_coefficient_a( aph2, -1, 3, l2) *  Jeta[l2-1])
                    / calc_coefficient_a( aph2, -1, 1, l2);
            for (int l2 = 0; l2 <= M-l1; ++l2){
                int aph3 = 2 * l1 + 2 * l2 - 1;
                Jzeta[1]  = calc_generalized_jacobi_polynomial( aph3, -1, 1, zeta);
                for (int l3 = 1; l3 < M-l1-l2; ++l3)
                    Jzeta[l3+1]  = ((zeta - calc_coefficient_a( aph3, -1, 2, l3)) *  Jzeta[l3] - calc_coefficient_a( aph3, -1, 3, l3) *  Jzeta[l3-1])
                        / calc_coefficient_a( aph3, -1, 1, l3);
                for (int l3 = 0; l3 <= M-l1-l2; ++l3){
                    Multiindex index_now = Unitary_Multiindex[0] * l1 + Unitary_Multiindex[1] * l2 + Unitary_Multiindex[2] * l3;
                    int ind_index = correspondence.index2number(index_now);
                    sol_[ind_index] = pow(1-y-z, l1) * Jxi[l1] * pow(1-z, l2) * Jeta[l2] * Jzeta[l3];
                }
            }
        }
        std::vector<double> sol_actual(n_index, 0);
        for (int ind_index = 0; ind_index < n_index; ++ind_index)
            for (int j = 0; j < n_transform_local[ind_index]; ++j)
                sol_actual[transform_local[ind_index][j]] += weight_transform_local[ind_index][j] * sol_[ind_index];
        for (int ind_index = 0; ind_index < n_index; ++ind_index)
            sol_sem_qpoint += sol_sem_local[ind_index] * sol_actual[ind_index];
        
        // update local error with linfinity norm
        double err_tmp = fabs(sol_qpoint - sol_sem_qpoint);
        err_local_l2 += Weight[i] * pow(err_tmp, 2);
    }
    return err_local_l2 * volume_local;
}
