#ifndef PhantomAnimator_class
#define PhantomAnimator_class

#include <iostream>
#include <ctime>
#include <functions.h>
#include <functional>
#include <iomanip>

#include <igl/readTGF.h>
#include <igl/writeTGF.h>
#include <igl/readDMAT.h>
#include <igl/writeDMAT.h>
#include <igl/writeMESH.h>
#include <igl/readMESH.h>
#include <igl/readPLY.h>
#include <igl/directed_edge_parents.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/boundary_conditions.h>
#include <igl/directed_edge_parents.h>
#include <igl/directed_edge_orientations.h>
#include <igl/deform_skeleton.h>
#include <igl/forward_kinematics.h>
#include <igl/doublearea.h>
#include <igl/dqs.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/normalize_row_sums.h>
#include <igl/Timer.h>
#include <igl/mat_max.h>

#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <Eigen/SparseCore>

#include "G4GeometryTolerance.hh"
#include "G4ThreeVector.hh"

class PhantomAnimator
{
//functions
public:
    PhantomAnimator();
    PhantomAnimator(string prefix, string profile);
    ~PhantomAnimator();

    void ReadTetMesh(string prefix);
    void PreparePhantom(string prefix);
    // bool ReadFiles(string prefix);
    // string CalibrateTo(string name);
    void Animate(RotationList vQ, Vector3d root);

    MatrixXd GetU(){ return U; }
    MatrixXi GetT(){ return T; }
    G4int GetNumOfTet(){return T.rows();}
    void GetTetVertices(G4int i, G4ThreeVector &a, G4ThreeVector &b, G4ThreeVector &c, G4ThreeVector &d)
    {
        Vector3d aV = U.row(T(i, 0));
        a = G4ThreeVector(aV(0), aV(1), aV(2));
        Vector3d bV = U.row(T(i, 1));
        b = G4ThreeVector(bV(0), bV(1), bV(2));
        Vector3d cV = U.row(T(i, 2));
        c = G4ThreeVector(cV(0), cV(1), cV(2));
        Vector3d dV = U.row(T(i, 3));
        d = G4ThreeVector(dV(0), dV(1), dV(2));
    }
    G4int GetMaterialIdx(G4int i) {return G4int(T(i, 4));}
    bool ReadProfileData(string profile);
    MatrixXd GetBox(){
        MatrixXd box(2, 3);
        box.row(0) = U.colwise().maxCoeff();
        box.row(1) = U.colwise().minCoeff();
        return box;
    }

//variables
private:
    MatrixXd C, V, U, W, Wj;
    MatrixXi BE, T;
    VectorXi P;
    // MatrixXd V_calib, C_calib;
    // vector<int> eyeIDs;
    vector<map<int, double>> cleanWeights;
    // map<int, double> lengths;

    // map<string, int> profileIDs;
    // vector<map<int, double>> jointLengths;
    // vector<Vector3d> eyeR_vec, eyeL_vec;

    // vector<int> CheckDegeneracy(const MatrixXd& VV, const MatrixXi& TT)
    // {
    //     double tol = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()*8./3.;
    
    //     MatrixXd A;
    //     igl::face_areas(VV, TT, A);
    //     VectorXd Amax = A.rowwise().maxCoeff();
    //     VectorXd vol;
    //     igl::volume(VV, TT, vol);
    //     vol = vol.array().abs();

    //     vector<int> degen;
    //     for(int i=0;i<T.rows();i++)
    //     {
    //         if(vol(i)<Amax(i)*tol)
    //             degen.push_back(i);
    //     }
    //     return degen;
    // }

};

#endif
