#include "PhantomAnimator.hh"
#include <igl/winding_number.h>
#include <igl/signed_distance.h>
// #include <igl/opengl/glfw/Viewer.h>
#include <igl/in_element.h>
#include <algorithm>

// PhantomAnimator::PhantomAnimator(string prefix)
//     : defaultPhantom("./phantoms/AM")
// {
//     if(!LoadPhantomWithWeightFiles(defaultPhantom)) LoadPhantom(defaultPhantom);
//     ReadProfileData("profile.txt");
// }
// PhantomAnimator &PhantomAnimator::Instance()
// {
//     static PhantomAnimator PhantomAnimator;
//     return PhantomAnimator;
// }

bool PhantomAnimator::LoadPhantom(string _phantom, bool isMale)
{  
    //try to load phantom with weight files
    //if it fails, generate weight files
    if(LoadPhantomWithWeightFiles(_phantom, isMale)) return true;
    
    // read phantom files
    string tgfFile;
    if(isMale) tgfFile = G4String(getenv("DCIR_PHANTOM_DIR")) + "AM.tgf";
    else tgfFile = G4String(getenv("DCIR_PHANTOM_DIR")) + "AF.tgf";
    cout << "Read " + tgfFile<< endl;
    if (!igl::readTGF(tgfFile, C, BE))
    {
        cout << "There is no " + _phantom + ".tgf";
        return false;
    }
    cout << "Read " + _phantom + ".ply" << endl;
    MatrixXd V_ply;
    MatrixXi F_ply;
    if (!igl::readPLY(_phantom + ".ply", V_ply, F_ply))
    {
        cout << "There is no " + _phantom + ".ply";
        return false;
    }

    // read apron files
    // cout << "Read " + _phantom + "_apron.ply" << endl;
    // igl::readPLY(_phantom + "_apron.ply", V_apron, F_apron);
    // U_apron = V_apron;

    // configure vertices inside the apron polygon
    // VectorXd S;
    // VectorXi I;
    // MatrixXd C0, N_tmp;
    // igl::signed_distance(V, V_apron, F_apron, igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL, S, I, C0, N_tmp);
    // vector<int> outside, inside;
    // for (int i = 0; i < S.rows(); i++)
    //     if (S(i) > 0)
    //         outside.push_back(i);
    //     else inside.push_back(i);
    // outApron.resize(outside.size());
    // outApron = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(outside.data(), outside.size());
    // apronMask = ArrayXd::Zero(V.rows());
    // VectorXd ones = VectorXd::Ones(outApron.rows());
    // igl::slice_into(ones, outApron, apronMask); // 0: inside apron, 1: outside apron

    // calculate weights
    cout << "generate bone points"<<endl;
    MatrixXd V1(V_ply.rows()+C.rows(), 3);
    V1<<V_ply, C;
    for(int i=0;i<BE.rows();i++)
    {
        RowVector3d vec = C.row(BE(i, 1))-C.row(BE(i, 0));
        int num = floor(vec.norm());
        V1.conservativeResize(V1.rows()+num-1, 3);
        RowVector3d start = C.row(BE(i, 0)) + vec/num;
        RowVector3d end = C.row(BE(i, 1)) - vec/num;
        V1.block(V1.rows()-num+1, 0, num-1, 1) = igl::LinSpaced<VectorXd>(num-1,start(0), end(0));
        V1.block(V1.rows()-num+1, 1, num-1, 1) = igl::LinSpaced<VectorXd>(num-1,start(1), end(1));
        V1.block(V1.rows()-num+1, 2, num-1, 1) = igl::LinSpaced<VectorXd>(num-1,start(2), end(2));
    }

    MatrixXd VT, WT_j, WT;
    MatrixXi TT, FT;
    cout << "perform tetrahedralization" << endl;
    igl::copyleft::tetgen::tetrahedralize(V1, F_ply, "pYqQ", VT, TT, FT);
    
    cout << "calculate joint BBW" << endl;
    MatrixXd bc = MatrixXd::Identity(C.rows(), C.rows());
    VectorXi b = igl::LinSpaced<VectorXi>(C.rows(), V.rows(), V.rows()+C.rows()-1);
    igl::BBWData bbw_data;
    bbw_data.active_set_params.max_iter = 10;
    bbw_data.verbosity = 0;
    if (!igl::bbw(VT, TT, b, bc, bbw_data, WT_j))
    {cout<<"BBW failed!"<<endl; return false;}    
    
    cout << "calculate bone BBW" << endl;
    MatrixXd bc1;
    igl::boundary_conditions(VT, TT, C, VectorXi(), BE, MatrixXi(), b, bc1);
    bc.resize(b.rows(), 18);
    bc = MatrixXd::Zero(b.rows(), 18);
    for(int i=0;i<BE.rows();i++)
    {
        bc.col(BE(i, 0)) += bc1.col(i);
    }
    if (!igl::bbw(VT, TT, b, bc, bbw_data, WT))
    {cout<<"BBW failed!"<<endl; return false;}
    // igl::writeDMAT(_phantom+".W0", WT, false);

    ReadNodeFile(_phantom + ".node");
    ReadEleFile(_phantom + ".ele");

    // gen bary.
    Eigen::SparseMatrix<double> bary;
    cout << "AABB->" << flush;
    igl::AABB<Eigen::MatrixXd, 3> tree;
    tree.init(VT, TT);
    Eigen::VectorXi I;
    cout << "in_element->" << flush;
    igl::in_element(VT, TT, V, tree, I);
    cout << "calculate vol." << endl;
    Eigen::VectorXd invol;
    igl::volume(VT, TT, invol);
    invol = -(invol * 6.).cwiseInverse();
    typedef Eigen::Triplet<double> T;
    std::vector<T> coeff;
    bary.resize(V.rows(), VT.rows());
    for (int i = 0; i < I.rows(); i++)
    {
      vector<Eigen::Vector3d> tet;
      for (int n = 0; n < 4; n++)
        tet.push_back(Eigen::Vector3d(VT.row(TT(I(i), n))));
      Eigen::Vector3d v = V.row(i);
      coeff.push_back(T(i, TT(I(i), 0), (v - tet[3]).cross(tet[1] - tet[3]).dot(tet[2] - tet[3]) * invol(I(i))));
      coeff.push_back(T(i, TT(I(i), 1), (tet[0] - tet[3]).cross(v - tet[3]).dot(tet[2] - tet[3]) * invol(I(i))));
      coeff.push_back(T(i, TT(I(i), 2), (tet[0] - tet[3]).cross(tet[1] - tet[3]).dot(v - tet[3]) * invol(I(i))));
      coeff.push_back(T(i, TT(I(i), 3), (tet[0] - v).cross(tet[1] - v).dot(tet[2] - v) * invol(I(i))));
      cout << "\rGenerating barycentric_coordinate..." << i << "/" << I.rows() - 1 << flush;
    }
    cout << endl;    
    bary.setFromTriplets(coeff.begin(), coeff.end());
    
    // get weight for full phantom
    Wj = bary * WT_j;
    WT = bary * WT;
    double epsilon = 0.01;
    Wj = Wj.array() * (Wj.array()>epsilon).cast<double>() ;
    igl::normalize_row_sums(Wj, Wj);
    WT = WT.array() * (WT.array()>epsilon).cast<double>() ;
    igl::normalize_row_sums(WT, WT);

    igl::writeDMAT(_phantom + ".Wj", Wj);
    igl::writeDMAT(_phantom + ".W", WT);

    //apron tetrahedralization
    // cout << "tetrahedralize apron" << endl;
    // VectorXi inApron(inside.size());
    // inApron = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(inside.data(), inside.size());
    // MatrixXd V_apron1(V_apron.rows()+inApron.rows(), 3);
    // V_apron1<<V_apron, igl::slice(V, inApron, 1); 
    // MatrixXd VT;
    // igl::copyleft::tetgen::tetrahedralize(V_apron1, F_apron, "pYqQ", VT, TT, FT);

    //apron bone weight
    // cout << "calculate bone BBW" << endl;
    // b = igl::LinSpaced<VectorXi>(inApron.rows(), V_apron.rows(), V_apron1.rows()-1);
    // bc = igl::slice(WT, inApron, 1);
    // double effW(0.1);
    // VectorXi effCol((bc.colwise().maxCoeff().array()>effW).cast<int>().sum());
    // for(int i=0, n=0;i<bc.cols();i++)
    //     if(bc.col(i).maxCoeff()>effW) effCol(n++) = i;
    // bc = igl::slice(bc, effCol, 2);
    // MatrixXd W_apron, W;
    // if (!igl::bbw(VT, TT, b, bc, bbw_data, W))
    // {cout<<"BBW failed!"<<endl; return false;}
    // W_apron = MatrixXd::Zero(V_apron.rows(), WT.cols());
    // igl::slice_into(W.topRows(V_apron.rows()), effCol, 2, W_apron);
    // W_apron = W_apron.array() * (W_apron.array()>epsilon).cast<double>() ;
    // igl::normalize_row_sums(W_apron, W_apron);
    // igl::writeDMAT(_phantom+".W_apron", W_apron, false);

    //apron bone weight
    // cout << "calculate joint BBW" << endl;
    // bc = igl::slice(Wj0, inApron, 1);
    // effCol.resize((bc.colwise().maxCoeff().array()>effW).cast<int>().sum());
    // for(int i=0, n=0;i<bc.cols();i++)
    //     if(bc.col(i).maxCoeff()>effW) effCol(n++) = i;
    // bc = igl::slice(bc, effCol, 2);
    // if (!igl::bbw(VT, TT, b, bc, bbw_data, W))
    // {cout<<"BBW failed!"<<endl; return false;}
    // Wj_apron = MatrixXd::Zero(V_apron.rows(), Wj0.cols());
    // igl::slice_into(W.topRows(V_apron.rows()), effCol, 2, Wj_apron);
    // Wj_apron = Wj_apron.array() * (Wj_apron.array()>epsilon).cast<double>() ;
    // igl::normalize_row_sums(Wj_apron, Wj_apron);
    // igl::writeDMAT(_phantom+".Wj_apron", Wj_apron, false);

    cleanWeights.clear();
    for (int i = 0; i < VT.rows(); i++)
    {
        double sum(0);
        map<int, double> vertexWeight;
        for (int j = 0; j < WT.cols(); j++)
        {
            if (WT(i, j) < epsilon)
                continue;
            vertexWeight[j] = WT(i, j);
            sum += WT(i, j);
        }
        for (auto &iter : vertexWeight)
            iter.second /= sum;
        cleanWeights.push_back(vertexWeight);
    }

    // cleanWeightsApron.clear();
    // for (int i = 0; i < W_apron.rows(); i++)
    // {
    //     double sum(0);
    //     map<int, double> vertexWeight;
    //     for (int j = 0; j < W_apron.cols(); j++)
    //     {
    //         if (W_apron(i, j) < epsilon)
    //             continue;
    //         vertexWeight[j] = W_apron(i, j);
    //         sum += W_apron(i, j);
    //     }
    //     for (auto &iter : vertexWeight)
    //         iter.second /= sum;
    //     cleanWeightsApron.push_back(vertexWeight);
    // }

    Initialize();
    return true;
}

bool PhantomAnimator::LoadPhantomWithWeightFiles(string _phantom, bool isMale)
{
    // read phantom files
    string tgfFile;
    if(isMale) tgfFile = G4String(getenv("DCIR_PHANTOM_DIR")) + "AM.tgf";
    else tgfFile = G4String(getenv("DCIR_PHANTOM_DIR")) + "AF.tgf";
    cout << "Read " + tgfFile << endl;
    if (!igl::readTGF(tgfFile, C, BE))
    {
        cout << "There is no " + _phantom + ".tgf";
        return false;
    }
    cout << "Read " + _phantom + ".W/Wj" << endl;
    MatrixXd W;
    if(!igl::readDMAT(_phantom+".W", W)) return false;
    if(!igl::readDMAT(_phantom+".Wj", Wj)) return false;

    ReadNodeFile(_phantom+".node");
    ReadEleFile(_phantom+".ele");

    // read apron files
    // cout << "Read " + _phantom + "_apron.ply" << endl;
    // igl::readPLY(_phantom + "_apron.ply", V_apron, F_apron);
    // U_apron = V_apron;

    // configure vertices inside the apron polygon
    // VectorXd S;
    // VectorXi I;
    // MatrixXd C0, N_tmp;
    // igl::signed_distance(V, V_apron, F_apron, igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL, S, I, C0, N_tmp);
    // vector<int> outside;
    // for (int i = 0; i < S.rows(); i++)
    //     if (S(i) > 0)
    //         outside.push_back(i);
    // outApron.resize(outside.size());
    // outApron = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(outside.data(), outside.size());
    // apronMask = ArrayXd::Zero(V.rows());
    // VectorXd ones = VectorXd::Ones(outApron.rows());
    // igl::slice_into(ones, outApron, apronMask); // 0: inside apron, 1: outside apron

    // cout<<"read "+_phantom+".W_apron"<<endl;
    // MatrixXd W_apron;
    // if(!igl::readDMAT(_phantom+".W_apron", W_apron))
    // {cout<<"There is no "+_phantom+".W_apron"; return false;}

    // cout<<"read "+_phantom+".Wj_apron"<<endl;
    // if(!igl::readDMAT(_phantom+".Wj_apron", Wj_apron))
    // {cout<<"There is no "+_phantom+".Wj_apron"; return false;}

    double epsilon(0.01);
    cleanWeights.clear();
    for (int i = 0; i < V.rows(); i++)
    {
        double sum(0);
        map<int, double> vertexWeight;
        for (int j = 0; j < W.cols(); j++)
        {
            if (W(i, j) < epsilon)
                continue;
            vertexWeight[j] = W(i, j);
            sum += W(i, j);
        }
        for (auto &iter : vertexWeight)
            iter.second /= sum;
        cleanWeights.push_back(vertexWeight);
    }

    // cleanWeightsApron.clear();
    // for (int i = 0; i < W_apron.rows(); i++)
    // {
    //     double sum(0);
    //     map<int, double> vertexWeight;
    //     for (int j = 0; j < W_apron.cols(); j++)
    //     {
    //         if (W_apron(i, j) < epsilon)
    //             continue;
    //         vertexWeight[j] = W_apron(i, j);
    //         sum += W_apron(i, j);
    //     }
    //     for (auto &iter : vertexWeight)
    //         iter.second /= sum;
    //     cleanWeightsApron.push_back(vertexWeight);
    // }

    Initialize();

    return true;
}
bool PhantomAnimator::Initialize()
{
    endC.clear();
    endC.resize(18);
    for(int i=0;i<BE.rows();i++)
     endC[BE(i, 0)].push_back(BE(i, 1));
    // additional variables
    // if(BE.rows()==0) igl::readTGF("./phantoms/AM.tgf", C, BE);
    igl::directed_edge_parents(BE, P);

    // distance to parent joint
    for (int i = 0; i < BE.rows(); i++)
        lengths[i] = (C.row(BE(i, 0)) - C.row(BE(i, 1))).norm();

    // joint number conv.
    vector<int> i2k = {0, 1, 2, 4, 5, 6, 7, 26, 11, 12, 13, 14, 18, 19, 20, 22, 23, 24};

    map<int, int> k2i;
    for (size_t i = 0; i < i2k.size(); i++)
        k2i[i2k[i]] = i;

    // preprocessing for KINECT data
    map<int, Quaterniond> kinectDataRot;
    vector<int> groups;
    Matrix3d tf2;
    tf2 << 0, -1, 0, 0, 0, -1, 1, 0, 0;
    for (int id : groups = {0, 1, 2, 18, 19, 20, 26})
        kinectDataRot[id] = Quaterniond(tf2);
    tf2 << 0, 1, 0, 0, 0, 1, 1, 0, 0;
    for (int id : groups = {22, 23, 24})
        kinectDataRot[id] = Quaterniond(tf2);
    tf2 << 0, 1, 0, 0, 0, -1, -1, 0, 0;
    for (int id : groups = {5, 6})
        kinectDataRot[id] = Quaterniond(tf2);
    tf2 << 0, -1, 0, 0, 0, 1, -1, 0, 0;
    for (int id : groups = {12, 13})
        kinectDataRot[id] = Quaterniond(tf2);
    tf2 << 1, 0, 0, 0, 0, 1, 0, -1, 0;
    kinectDataRot[11] = Quaterniond(tf2);
    tf2 << 1, 0, 0, 0, 0, -1, 0, 1, 0;
    kinectDataRot[4] = Quaterniond(tf2);
    tf2 << 0, 1, 0, 1, 0, 0, 0, 0, -1;
    kinectDataRot[7] = Quaterniond(tf2);
    tf2 << 0, -1, 0, 1, 0, 0, 0, 0, 1;
    kinectDataRot[14] = Quaterniond(tf2);

    map<int, Vector3d> desiredOrt;
    groups = {12, 13, 5, 6, 22, 23, 18, 19};
    for (int id : groups)
        desiredOrt[id] = Vector3d(0, 1, 0);
    desiredOrt[4] = Vector3d(1, 0, 0);
    desiredOrt[11] = Vector3d(-1, 0, 0);
    // for (int i = 0; i < BE.rows(); i++)
    // {
    //     if (desiredOrt.find(i2k[BE(i, 0)]) == desiredOrt.end())
    //         alignRot.push_back(kinectDataRot[i2k[BE(i, 0)]]);
    //     else
    //     {
    //         Vector3d v = (C.row(k2i[i2k[BE(i, 0)] + 1]) - C.row(BE(i, 0))).transpose();
    //         alignRot.push_back(kinectDataRot[i2k[BE(i, 0)]] * GetRotMatrix(v, desiredOrt[i2k[BE(i, 0)]]));
    //     }
    // }

    alignRot.resize(18);
    for(int i=0;i<alignRot.size();i++)
    {
        if (desiredOrt.find(i2k[i]) == desiredOrt.end())
            alignRot[i] = kinectDataRot[i2k[i]];
        else
        {
            /////careful!!
            Vector3d v = (C.row(i + 1) - C.row(i)).transpose();
            alignRot[i] = kinectDataRot[i2k[i]] * GetRotMatrix(v, desiredOrt[i2k[i]]);
        }
    }

    // cout << "generating F_area, F_incident, F_area_eye matrix..." << flush;
    // VectorXd F_area;
    // igl::doublearea(V, F, F_area);
    // F_area.cwiseSqrt();
    // VectorXd W_avgSkin = ArrayXd::Zero(V.rows());

    // for (int i = 0; i < F.rows(); i++)
    // {
    //     for (int j = 0; j < 3; j++)
    //         W_avgSkin(F(i, j)) += F_area(i);
    // }
    // W_avgSkin = W_avgSkin / W_avgSkin.sum();

    // VectorXd W_Lens(eye2ply.size());
    // for (size_t i = 0; i < eye2ply.size(); i++)
    //     W_Lens(i) = W_avgSkin(eye2ply[i]);

    // W_avgSkin = W_avgSkin.array() / W_avgSkin.sum();
    // W_Lens = W_Lens.array() / W_Lens.sum();

    cout << "done" << endl;
    // V_calib = V;
    // C_calib = C;
    // V_calib_apron = V_apron;
    return true;
}

bool PhantomAnimator::CalibrateTo(string name)
{
    if(profileIDs.find(name)==profileIDs.end()) return false;
    int id = profileIDs[name];
    map<int, double> calibLengths = jointLengths[id];
    RowVector3d eyeL_pos = eyeL_vec[id];
    RowVector3d eyeR_pos = eyeR_vec[id];

    MatrixXd jointTrans = MatrixXd::Zero(C.rows(), 3);
    int headJ(24), eyeLJ(22), eyeRJ(23);
    for (int i = 0; i < BE.rows(); i++)
    {
        if (calibLengths.find(i) == calibLengths.end())
        {
            calibLengths[i] = lengths[i];
            jointTrans.row(BE(i, 1)) = jointTrans.row(BE(P(i), 1));
            continue;
        }
        double ratio = calibLengths[i] / lengths[i];
        cout << i << " : " << lengths[i] << " -> " << calibLengths[i] << " (" << ratio * 100 << " %)" << endl;
        jointTrans.row(BE(i, 1)) = (1 - ratio) * (C.row(BE(i, 0)) - C.row(BE(i, 1)));
        if (P(i) < 0)
            continue;
        jointTrans.row(BE(i, 1)) += jointTrans.row(BE(P(i), 1));
    }

    jointTrans.row(eyeLJ) = C.row(headJ) + jointTrans.row(headJ) + eyeL_pos - C.row(eyeLJ);
    jointTrans.row(eyeRJ) = C.row(headJ) + jointTrans.row(headJ) + eyeR_pos - C.row(eyeRJ);
    jointTrans.row(headJ) = MatrixXd::Zero(1, 3);
    jointTrans(headJ, 1) = (jointTrans(eyeLJ, 1) + jointTrans(eyeRJ, 1)) * 0.5;
    C = C + jointTrans;

    cout << Wj.rows() << "*" << Wj.cols() << endl;
    cout << jointTrans.rows() << "*" << jointTrans.cols() << endl;
    V = V + Wj * jointTrans.block(0, 0, C.rows() - 1, 3);

    // cout << V_apron.rows() << " " << V_apron.cols() << endl;
    // cout << Wj_apron.rows() << " " << Wj_apron.cols() << endl;
    // cout << C.rows() - 1 << endl;

    // V_apron = V_apron + Wj_apron * jointTrans.block(0, 0, C.rows() - 1, 3);
    // MatrixXd jt = jointTrans.block(0,0,C.rows()-1,3);
    return true;
}

void PhantomAnimator::Animate(RotationList vQ, RowVector3d rootC)
{
    vector<Vector3d> vT;
    MatrixXd C_new = C;
    C_new.row(0) = rootC; // set root
    for (int i = 0; i < 18;i++)//BE.rows(); i++)
    {
        vQ[i] = vQ[i] * alignRot[i];
        Affine3d a;
        a = Translation3d(Vector3d(C_new.row(i).transpose())) * vQ[i].matrix() * Translation3d(Vector3d(-C.row(i).transpose()));
        vT.push_back(a.translation());
        for(int e:endC[i]) C_new.row(e) = a * Vector3d(C_new.row(e));
    }
    myDqs(V, cleanWeights, vQ, vT, U);
    myDqs(V_apron, cleanWeightsApron, vQ, vT, U_apron);
}

bool PhantomAnimator::ReadProfileData(string fileName)
{
    ifstream ifs(fileName);
    if (!ifs.is_open())
    {
        cout << "There is no " + fileName << endl;
        return false;
    }
    profileIDs.clear();
    jointLengths.clear();
    eyeR_vec.clear();
    eyeL_vec.clear();

    int num;
    ifs >> num;
    string firstLine;
    getline(ifs, firstLine);
    jointLengths.resize(num);
    eyeR_vec.resize(num);
    eyeL_vec.resize(num);

    getline(ifs, firstLine);
    stringstream ss(firstLine);
    vector<int> boneIDs;
    for (int i = 0; i < 17; i++)
    {
        int boneID;
        ss >> boneID;
        boneIDs.push_back(boneID);
    }

    for (int i = 0; i < num; i++)
    {
        string aLine;
        getline(ifs, aLine);
        stringstream ss1(aLine);
        string name;
        double d;
        ss1 >> name;
        profileIDs[name] = i;
        for (int j = 0; j < 17; j++)
        {
            double d;
            ss1 >> d;
            jointLengths[i][boneIDs[j]] = d;
        }
        for (int j = 0; j < 3; j++)
        {
            ss1 >> d;
            eyeL_vec[i](j) = d;
        }
        for (int j = 0; j < 3; j++)
        {
            ss1 >> d;
            eyeR_vec[i](j) = d;
        }
    }

    ifs.close();
    return true;
}

void PhantomAnimator::ReadNodeFile(G4String nodefile)
{
    ifstream ifsN(nodefile);
    if (!ifsN.is_open())
    {
        G4Exception("PhantomAnimator::ReadNodeFile","",FatalErrorInArgument,
	          G4String("There is no " + nodefile).c_str());
    }
    cout << "reading node/ele file" << endl;
    int num;
    G4String firstlineN;
    getline(ifsN, firstlineN);
    stringstream ss(firstlineN);
    ss >> num;
    V.resize(num, 3);
    for (int i = 0; i < num; i++)
    {
      int tmp;
      ifsN >> tmp >> V(i, 0) >> V(i, 1) >> V(i, 2);
    }
    U=V;
    ifsN.close();
}

void PhantomAnimator::ReadEleFile(G4String eleFile)
{
    ifstream ifsE(eleFile);
     if (!ifsE.is_open())
     {
                G4Exception("PhantomAnimator::ReadEleFile","",FatalErrorInArgument,
	          G4String("There is no " + eleFile).c_str());
     }
     string firstlineE;
     getline(ifsE, firstlineE);
     stringstream ss1(firstlineE);
     int num;
     ss1 >> num;
     T.resize(num, 5);
     for (int i = 0; i < num; i++)
     {
       int tmp;
       ifsE >> tmp >> T(i, 0) >> T(i, 1) >> T(i, 2) >> T(i, 3) >> T(i, 4);
     }
     ifsE.close();    
}