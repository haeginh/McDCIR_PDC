#ifndef WEBServerConnect_HH
#define WEBServerConnect_HH 1

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cstring>
#include <iomanip>

#include <mysql.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <sys/ioctl.h>
#include <netinet/in.h>
#include <net/if.h>

#include "functions.h"
#include "TETModelImport.hh"

using namespace std;
using namespace Eigen;

enum QUERY { SUCCESS, FAIL };
enum CALC_STATUS { NOTinProgress, InProgress, Done };
enum SERV_STATUS { ALLOCATED, EMPTY };
typedef pair<string, int> COLIDX;

class WEBServerConnect
{
public:
    WEBServerConnect(int id);
    ~WEBServerConnect();

    // Processing web-server
    int ReadCalculationNumber();
    SERV_STATUS ProcessData();
    void UpdateFlagToDone();
    void ResetData();

    // Access Functions
    // Execution data
    string GetMyIPAddress();
    int GetThreadNo() { return threadNo; }
    int GetNPS()      { return NPS;      }
    MYSQL* GetConnection() { return conn; }
    vector<COLIDX> GetColDoseResults() { return col_doseResults; }
    int GetCalculationFrameNo() { return calNo; }

    // Tracking data
    string GetHospitalName() { return hospital; }
    int GetFrameNo() { return frameNo; }
    double GetFrameTime() { return frameTime; }
    int GetTubeVoltage() { return tubeVoltage; }
    int GetTubeCurrent() { return tubeCurrent; }
    double GetDAP() { return DAP; }
    int GetFD() { return FD; }
    int GetSID() { return SID; }
    int GetRAOLAO() { return RAOLAO; }
    int GetCRANCAUD() { return CRANCAUD; }
    Vector3d GetTableTranslation() { return table_trans; }
    double GetTableTheta() { return table_theta; }
    Vector3d GetGlassTranslation() { return glass_trans; }
    Quaterniond GetGlassQuaternion() { return glass_quat; }
    Vector3d GetDoctorRoot() { return doctor_root; }
    RotationList GetDoctorRotationList() { return doctor_vQ; }

    // Dose data
    void SendDoseResultsToWebServer(map<G4int, G4String> nameMap,
                                    map<G4int, pair<G4double, G4double>> doses,
                                    pair<G4double, G4double> effective_DRF);

private:
    void ReadDBTableColumn(string fileName);
    void ReadProgramParameter(string fileName, int id);
    void DecryptFile(string fileName);
    void finish_with_error(MYSQL* conn) {
        cerr << "SQL Error !! - " + string(mysql_error(conn)) << endl; 
        mysql_close(conn); 
        exit(1);
    }
    


private:
    bool isFirst;

    string DB_HOST, DB_USER, DB_PASS, DB_NAME, DB_TABLE, MyIP;
    int PORT_ID;

    MYSQL* conn;
    MYSQL_RES* result;
    MYSQL_ROW row;

    // indices
    COLIDX col_num, col_frameNo, col_flag, col_server, 
           col_hospital, col_time, col_frameTime,
           col_tVolt, col_tCurr, col_DAP, col_FD, col_SID, col_RAOLAO, col_CRANCAUD;
    vector<COLIDX> col_table_xyztheta, col_glass_xyz, col_glass_wxyz, 
                   col_radiologist_xyz, col_radiologist_wxyz, col_doseResults;
           
    
    // Execution data
    int threadNo, NPS;
    int calNo;
    // Tracking data
    // 1. Frame
    string hospital;
    int frameNo;
    double frameTime;
    // 2. C-arm
    int tubeVoltage;
    int tubeCurrent;
    double DAP;
    int FD;
    int SID;
    int RAOLAO;
    int CRANCAUD;
    // 3. Operating table
    Vector3d table_trans;
    double table_theta;
    // 4. Pb glass
    Vector3d glass_trans;
    Quaterniond glass_quat;
    // 5. Doctor
    Vector3d doctor_root;
    RotationList doctor_vQ;
        
};


#endif