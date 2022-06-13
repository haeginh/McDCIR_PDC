#include "WEBServerConnect.hh"
#include "G4SystemOfUnits.hh"

WEBServerConnect::WEBServerConnect(int id): 
isFirst(true)
{
    string WebServerInfoFile = "PDC.dat";
    string DBTableColumnFile = "PDC.column";
    string ServerConfigFile  = "PDC.server";

    conn = mysql_init(NULL);
    DecryptFile(WebServerInfoFile);
    ReadDBTableColumn(DBTableColumnFile);
    ReadProgramParameter(ServerConfigFile, id);
    MyIP = GetMyIPAddress();
    
    // Connect SQL server
	if( mysql_real_connect(conn, DB_HOST.c_str(), DB_USER.c_str(), DB_PASS.c_str(), DB_NAME.c_str(), PORT_ID, NULL, 0) != NULL ) {
		cout << "=======================================================================" << endl;
		cout << "Web server is successfully connected!" << endl;
		cout << "DB Host IP: " << DB_HOST << endl;
		cout << "DB User   : " << DB_USER << endl;
		cout << "DB Name   : " << DB_NAME << endl;
		cout << "DB TABLE  : " << DB_TABLE << endl;
		cout << "Server IP : " << MyIP << endl;
		cout << "Port ID   : " << PORT_ID << endl;
		cout << "=======================================================================" << endl;
	} else finish_with_error(conn);
}

WEBServerConnect::~WEBServerConnect()
{
    cout << "Web server is successfully disconnected!" << endl;
    mysql_close(conn);
}

void WEBServerConnect::SendDoseResultsToWebServer(map<G4int, G4String> nameMap,
                                                  map<G4int, pair<G4double, G4double>> doses,
                                                  pair<G4double, G4double> effective_DRF)
{    
    vector<double> doseResults;
    if (isFirst) {
        int idx(0);
        for (auto itr: nameMap) {
            if (itr.first == -1 || itr.first == -2) continue;
            col_doseResults[idx++].first  = itr.second;
            col_doseResults[idx++].first  = itr.second + "_err";
            doseResults.push_back(doses[itr.first].first  );
            doseResults.push_back(doses[itr.first].second );
        }
        col_doseResults[idx+0].first = "HandL";
        col_doseResults[idx+1].first = "HandL_err";
        col_doseResults[idx+2].first = "HandR";
        col_doseResults[idx+3].first = "HandR_err";
        col_doseResults[idx+4].first = "Eff(DRF)";
        col_doseResults[idx+5].first = "Eff(DRF)_err";
        doseResults.push_back(0.);
        doseResults.push_back(0.);
        doseResults.push_back(0.);
        doseResults.push_back(0.);
        doseResults.push_back(effective_DRF.first);
        doseResults.push_back(effective_DRF.second);
        isFirst = false;
    }

    string sql;
    for (size_t i=0; i<col_doseResults.size(); i++) {        
        double result = doseResults[i];
        if ( i % 2 == 1 ) 
            result = result * (joule/kg) / 1e12;
        if (isnan(result)) result = 0.;
        sql = "UPDATE `" + DB_TABLE + "` SET `" + col_doseResults[i].first + "` = '" + to_string(result/(joule/kg)*1e12)
            + "' Where `" + col_num.first + "` = '" + to_string(calNo) + "';";
        cout << "SQL query => " << sql << endl;
        if(mysql_query(conn, sql.c_str()) != QUERY::SUCCESS) finish_with_error(conn);
    }
}

int WEBServerConnect::ReadCalculationNumber()
{
    string sql;
    sql = "SELECT `" + col_flag.first + "` FROM `" + DB_TABLE + "` WHERE `" + col_num.first + "` = -1;";
    if(mysql_query(conn, sql.c_str()) != QUERY::SUCCESS) finish_with_error(conn);
    result = mysql_store_result(conn);
    row = mysql_fetch_row(result);
    calNo = atoi(row[0]);
    return calNo;
}

SERV_STATUS WEBServerConnect::ProcessData()
{
    string sql;
    // Update total frame # flag
    sql = "UPDATE `" + DB_TABLE + "` SET `" + col_flag.first 
        + "` = '" + to_string(calNo + 1) + "' Where `" + col_num.first + "` = -1;";
    cout << "SQL query => " << sql << endl;
    if(mysql_query(conn, sql.c_str()) != QUERY::SUCCESS) finish_with_error(conn);

    // Access the frame data
    sql = "SELECT * FROM `" + DB_TABLE + "` WHERE Number = '" + to_string(calNo) + "';";
    cout << "SQL query => " << sql << endl;
    if(mysql_query(conn, sql.c_str()) != QUERY::SUCCESS) finish_with_error(conn);
    result = mysql_store_result(conn);
    row = mysql_fetch_row(result);

    // Check the execution flag status
    if (atoi(row[col_flag.second]) != (int)CALC_STATUS::NOTinProgress) {
        cerr << "WARNING !! - The frame is already allocated in other server." << endl;
        return SERV_STATUS::ALLOCATED;
    }

    // Update parameters
    //
    // ====== Execution ======
    // 1. Flag
    sql = "UPDATE `" + DB_TABLE + "` SET `" + col_flag.first 
        + "` = '" + to_string(CALC_STATUS::InProgress) + "' Where `" + col_num.first + "` = '" + to_string(calNo) + "';";
    cout << "SQL query => " << sql << endl;
    if(mysql_query(conn, sql.c_str()) != QUERY::SUCCESS) finish_with_error(conn);

    // 2. Server IP
    sql = "UPDATE `" + DB_TABLE + "` SET `" + col_server.first + "` = '" + MyIP 
        + "' Where `" + col_num.first + "` = '" + to_string(calNo) + "';";
    cout << "SQL query => " << sql << endl;
    if(mysql_query(conn, sql.c_str()) != QUERY::SUCCESS) finish_with_error(conn);

    // 3. Time
    time_t timer = time(NULL);
    struct tm* t = localtime(&timer);
    string timeStr = to_string(t->tm_year + 1900) + "." + to_string(t->tm_mon + 1) + "." + to_string(t->tm_mday) 
           + " / " + to_string(t->tm_hour) + ":" + to_string(t->tm_min) + ":" + to_string(t->tm_sec);
    sql = "UPDATE `" + DB_TABLE + "` SET `" + col_time.first + "` = '" + timeStr 
        + "' Where `" + col_num.first + "` = '" + to_string(calNo) + "';";
    cout << "SQL query => " << sql << endl;
    if(mysql_query(conn, sql.c_str()) != QUERY::SUCCESS) finish_with_error(conn);
    // ====== Execution ======

    // ====== Tracking ======
    // 1. Frame
    hospital  = string(row[col_hospital.second]);
    frameNo   = atof( row[col_frameNo.second] );
    frameTime = atof( row[col_frameTime.second] );

    // 2. C-arm
    tubeVoltage = atof( row[col_tVolt.second] );
    tubeCurrent = atof( row[col_tCurr.second] );
    DAP = atof( row[col_DAP.second] );
    FD  = atof( row[col_FD.second] );
    SID = atof( row[col_SID.second] );
    RAOLAO = atof( row[col_RAOLAO.second] );
    CRANCAUD = atof( row[col_CRANCAUD.second] );

    // 3. Operatin Table
    table_trans = Vector3d( atof( row[col_table_xyztheta[0].second]), 
                            atof( row[col_table_xyztheta[1].second]), 
                            atof( row[col_table_xyztheta[2].second]) );
    table_theta = atof( row[col_table_xyztheta[3].second]);

    // 4. Pb Glass
    glass_trans = Vector3d( atof( row[col_glass_xyz[0].second]), 
                            atof( row[col_glass_xyz[1].second]), 
                            atof( row[col_glass_xyz[2].second]) );
    glass_quat = Quaterniond( atof( row[col_glass_wxyz[0].second]), 
                              atof( row[col_glass_wxyz[1].second]), 
                              atof( row[col_glass_wxyz[2].second]), 
                              atof( row[col_glass_wxyz[3].second]) );

    // 5. Doctor
    doctor_root = Vector3d( atof(row[col_radiologist_xyz[0].second]),
                            atof(row[col_radiologist_xyz[1].second]),
                            atof(row[col_radiologist_xyz[2].second]) );
    for (int i=0; i<22; i++) {
        double w = atof( row[col_radiologist_wxyz[ 4 * i + 0 ].second] );
        double x = atof( row[col_radiologist_wxyz[ 4 * i + 1 ].second] );
        double y = atof( row[col_radiologist_wxyz[ 4 * i + 2 ].second] );
        double z = atof( row[col_radiologist_wxyz[ 4 * i + 3 ].second] );
        doctor_vQ.push_back(Quaterniond(w,x,y,z));
    }
    // ====== Tracking ======

    cout << "Reading the data of number '" << row[col_num.second] << "'..." << endl;
    return SERV_STATUS::EMPTY;
}

void WEBServerConnect::UpdateFlagToDone()
{
    string sql;
    // Update the flag
    sql = "UPDATE `" + DB_TABLE + "` SET `" + col_flag.first + "` = '" + to_string(CALC_STATUS::Done) 
        + "' Where `" + col_num.first + "` = '" + to_string(calNo) + "';";
    cout << "SQL query => " << sql << endl;
    if(mysql_query(conn, sql.c_str()) != QUERY::SUCCESS) finish_with_error(conn);


    // Clear the frmae data
    doctor_root = Vector3d::Zero();
    doctor_vQ.clear();
    frameTime = 0.;
    tubeVoltage = 0;
    tubeCurrent = 0;
    DAP = 0.;
    FD = 0;
    SID = 0;
    RAOLAO = 0;
    CRANCAUD = 0;
    table_trans = Vector3d::Zero();
    table_theta = 0;
    glass_trans = Vector3d::Zero();
    glass_quat = Quaterniond::Identity();
}



void WEBServerConnect::ReadDBTableColumn(string fileName)
{
    cout << "Read DB Table Column" << endl;
    ifstream ifs(fileName);
    if (!ifs.is_open()) { 
        cerr << fileName << " was not opened" << endl; exit(1);
    }
    string dump;
    while (getline(ifs, dump)) {
        stringstream ss(dump);
        ss >> dump;
        if      (dump == "Number")      { col_num.first       = dump; ss >> col_num.second;       } 
        else if (dump == "Flag")        { col_flag.first      = dump; ss >> col_flag.second;      } 
        else if (dump == "Hospital")    { col_hospital.first  = dump; ss >> col_hospital.second;  } 
        else if (dump == "Server")      { col_server.first    = dump; ss >> col_server.second;    } 
        else if (dump == "Time")        { col_time.first      = dump; ss >> col_time.second;      } 
        else if (dump == "FrameNo")     { col_frameNo.first   = dump; ss >> col_frameNo.second;   } 
        else if (dump == "FrameTime")   { col_frameTime.first = dump; ss >> col_frameTime.second; }
        else if (dump == "TubeVoltage") { col_tVolt.first     = dump; ss >> col_tVolt.second;     }
        else if (dump == "TubeCurrent") { col_tCurr.first     = dump; ss >> col_tCurr.second;     }
        else if (dump == "DAP")         { col_DAP.first       = dump; ss >> col_DAP.second;       }
        else if (dump == "FD")          { col_FD.first        = dump; ss >> col_FD.second;        }
        else if (dump == "SID")         { col_SID.first       = dump; ss >> col_SID.second;       }
        else if (dump == "RAO/LAO")     { col_RAOLAO.first    = dump; ss >> col_RAOLAO.second;    }
        else if (dump == "CRAN/CAUD")   { col_CRANCAUD.first  = dump; ss >> col_CRANCAUD.second;  }
        else if (dump == "Table_xyztheta") {
            for (int i=0; i<4; i++) {
                ss >> dump;
                string xyztheta;
                if      (i == 0) xyztheta = "Table_x";
                else if (i == 1) xyztheta = "Table_y";
                else if (i == 2) xyztheta = "Table_z";
                else if (i == 3) xyztheta = "Table_theta";
                col_table_xyztheta.push_back( make_pair(xyztheta, stoi(dump) ));
            }
        }
        else if (dump == "Glass_xyz") {
            for (int i=0; i<3; i++) {
                ss >> dump;
                string xyz;
                if      (i == 0) xyz = "Glass_tx";
                else if (i == 1) xyz = "Glass_ty";
                else if (i == 2) xyz = "Glass_tz";
                col_glass_xyz.push_back( make_pair(xyz, stoi(dump) ));
            }
        }
        else if (dump == "Glass_wxyz") {
            for (int i=0; i<4; i++) {
                ss >> dump;
                string wxyz;
                if      (i == 0) wxyz = "Glass_qw";
                else if (i == 1) wxyz = "Glass_qx";
                else if (i == 2) wxyz = "Glass_qy";
                else if (i == 3) wxyz = "Glass_qz";
                col_glass_wxyz.push_back( make_pair(wxyz, stoi(dump) ));
            }
        }
        else if (dump == "Radiologist_xyz") {
            for (int i=0; i<3; i++) {
                ss >> dump;
                string xyz;
                if      (i == 0) xyz = "Doc_x";
                else if (i == 1) xyz = "Doc_y";
                else if (i == 2) xyz = "Doc_z";
                col_radiologist_xyz.push_back( make_pair(xyz, stoi(dump) ));
            }
        }
        else if (dump == "Radiologist_wxyz") {
            ss >> dump;
            int n = dump.find("-");
            int start = stoi(dump.substr(0,n));
            int end   = stoi(dump.substr(n+1,dump.size()));
            for (int i=start; i<=end; i++) {
                int idx = (i-start) * 0.25;
                string wxyz;
                if      ( (i-start) % 4 == 0 ) wxyz = "Doc_w" + to_string(idx);
                else if ( (i-start) % 4 == 1 ) wxyz = "Doc_x" + to_string(idx);
                else if ( (i-start) % 4 == 2 ) wxyz = "Doc_y" + to_string(idx);
                else if ( (i-start) % 4 == 3 ) wxyz = "Doc_z" + to_string(idx);
                col_radiologist_wxyz.push_back( make_pair(wxyz, i) );
            }
                
        }
        else if (dump == "DoseResults") {
            ss >> dump;
            int n = dump.find("-");
            int start = stoi(dump.substr(0,n));
            int end   = stoi(dump.substr(n+1,dump.size()));
            for (int i=start; i<=end; i++) {
                col_doseResults.push_back( make_pair("Dose", i) );
            }
        }
    }
}

void WEBServerConnect::ReadProgramParameter(string fileName, int id)
{
    ifstream ifs(fileName);
    if (!ifs.is_open()) { 
        cerr << fileName << " was not opened" << endl; exit(1);
    }

    string dump;
    getline(ifs, dump);
    bool idChk(false);
    while (getline(ifs, dump)) {
        stringstream ss(dump);
        ss >> dump;
        if (stoi(dump) == id) {
            ss >> threadNo >> NPS;
            idChk = true;
            break;
        }
    }
    if (!idChk) { cerr << "Server ID '"<< id << "' was not configured" << endl; exit(1); }
}


void WEBServerConnect::DecryptFile(string fileName)
{
    char ch;
    fstream fps, fpt;
    fps.open("tmp.txt", fstream::out);
    if(!fps) {
        cout << "Error Occured, Opening the Source File" << endl;
        exit(1);
    }
    fpt.open(fileName, fstream::in);
    if(!fpt) {
        cout << "Error Occured while Opening/Creating tmp file!" << endl;
        exit(1);
    }
    while (fpt >> noskipws >> ch) {
        ch = ch - 100;
        fps << ch;
    }
    fps.close();
    fpt.close();
    cout << "File '" << fileName << "' is sucessfully decrtyped!" << endl;
    ifstream ifs("tmp.txt");
    string dump;
    ifs >> dump >> DB_HOST;
	ifs >> dump >> DB_USER;
	ifs >> dump >> DB_PASS;
	ifs >> dump >> DB_NAME;
	ifs >> dump >> DB_TABLE;
	ifs >> dump >> PORT_ID;
    ifs.close();
    remove("tmp.txt");

    return;
}

string WEBServerConnect::GetMyIPAddress()
{
	char myip[20];
	int sockfd = socket(AF_INET, SOCK_DGRAM, 0);

	const char* kGoogleDnsIp = "8.8.8.8";
	int kDnsPort = 53;

	struct sockaddr_in serv;
	struct sockaddr_in host_name;

	memset(&serv, 0, sizeof(serv));

	serv.sin_family = AF_INET;
	serv.sin_addr.s_addr = inet_addr(kGoogleDnsIp);
	serv.sin_port = htons(kDnsPort);

	if( connect(sockfd, (struct sockaddr *)&serv, sizeof(serv)) < 0 ) printf("[-] sock connect for get ipaddr faild!\n");

	socklen_t host_len = sizeof(host_name);
	if( getsockname(sockfd, (struct sockaddr *)&host_name, &host_len) < 0 ) printf("[-] getsockname faild!\n");

	inet_ntop(AF_INET, &host_name.sin_addr, myip, sizeof(myip));
	close(sockfd);

	return string(myip);
}

void WEBServerConnect::ResetData()
{
    string sql;
    sql = "SELECT * FROM " + DB_TABLE + ";";
    if (mysql_query(conn, sql.c_str()) != QUERY::SUCCESS) finish_with_error(conn);
    result = mysql_store_result(conn);
    while ( (row = mysql_fetch_row(result)) != NULL )
    {
        int calNo = atoi(row[0]);
        // Update program parameter
        //
        // 1. Flag
        sql = "UPDATE `" + DB_TABLE + "` SET `" + col_flag.first + "` = '" + to_string(CALC_STATUS::NOTinProgress) 
            + "' Where `" + col_num.first + "` = '" + to_string(calNo) + "';";
        cout << "SQL query => " << sql << endl;
        if(mysql_query(conn, sql.c_str()) != QUERY::SUCCESS) finish_with_error(conn);

        // 2. Hospital
        sql = "UPDATE `" + DB_TABLE + "` SET `" + col_hospital.first + "` = '" + "-" 
            + "' Where `" + col_num.first + "` = '" + to_string(calNo) + "';";
        cout << "SQL query => " << sql << endl;
        if(mysql_query(conn, sql.c_str()) != QUERY::SUCCESS) finish_with_error(conn);

        // 3. Server IP
        sql = "UPDATE `" + DB_TABLE + "` SET `" + col_server.first + "` = '" + "-" 
            + "' Where `" + col_num.first + "` = '" + to_string(calNo) + "';";
        cout << "SQL query => " << sql << endl;
        if(mysql_query(conn, sql.c_str()) != QUERY::SUCCESS) finish_with_error(conn);

        // 4. Time
        sql = "UPDATE `" + DB_TABLE + "` SET `" + col_time.first + "` = '" + "-" 
            + "' Where `" + col_num.first + "` = '" + to_string(calNo) + "';";
        cout << "SQL query => " << sql << endl;
        if(mysql_query(conn, sql.c_str()) != QUERY::SUCCESS) finish_with_error(conn);
    }
}