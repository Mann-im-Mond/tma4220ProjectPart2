#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>

using namespace std;

int main(int argc, char **argv) {
	if(argc < 2) {
		cerr << "File usage: " << argv[0] << " <inputfile> <outputfile>" << endl;
		return 1;
	}

	string inFileName(argv[1]);
	if(inFileName.substr(inFileName.length()-4).compare(".msh") != 0) {
		cerr << "Input file should be a .msh file" << endl;
		exit(2);
	}

	string outFileName;
	if(argc < 3) {
		outFileName = inFileName.substr(0,inFileName.length()-4);
	} else {
		outFileName = argv[2]; 
	}

	ifstream ifs(inFileName.c_str(), ifstream::in);
	if(!ifs.is_open()) {
		cerr << "Error opening file " << inFileName << endl;
		exit(3);
	}
	
	string nodeFileName = outFileName + "_nodes.m";
	string ptsFileName  = outFileName + "_pts.m";
	string lineFileName = outFileName + "_lines.m";
	string triFileName  = outFileName + "_tri.m";
	string tetrFileName = outFileName + "_tetr.m";
	ofstream nodeOut(nodeFileName.c_str(), ofstream::out);
	ofstream ptsOut(ptsFileName.c_str(),   ofstream::out);
	ofstream lineOut(lineFileName.c_str(), ofstream::out);
	ofstream triOut( triFileName.c_str(),  ofstream::out);
	ofstream tetrOut(tetrFileName.c_str(), ofstream::out);

	string   *allFileName[] = {&nodeFileName, &ptsFileName, &lineFileName, &triFileName, &tetrFileName};
	ofstream *allOut[]      = {&nodeOut,      &ptsOut,      &lineOut,      &triOut,      &tetrOut     };
	int nOutputFiles       = 5;

	for(int i=0; i<nOutputFiles; i++) {
		if(!allOut[i]->is_open() ) {
			cerr << "Error opening file " << *allFileName[i] << endl;
			exit(3);
		}
	}
	
	int n ;
	bool isDone = false;
	while(!isDone) {
		string inLine;
		getline(ifs, inLine);

		if(inLine.compare("$MeshFormat") == 0) {
			getline(ifs, inLine);
			if(inLine.compare("2.2 0 8") != 0) {
				cerr << "Unsupported msh format. Only version 2.2 supported" << endl;
				exit(4);
			}
			getline(ifs, inLine); // should be $EndMeshFormat, don't bother checking

		} else if(inLine.compare("$Nodes") == 0) {
			ifs >> n;
			getline(ifs, inLine); // just the endline character
			for(int i=0; i<n; i++) {
				getline(ifs, inLine);
				nodeOut << inLine << endl;
			}
			getline(ifs, inLine); // should be $EndNodes, but don't bother checking
		} else if(inLine.compare("$Elements") == 0) {
			int elmID;
			int nParam;
			int elmType;
			vector<int> metaInfo;
			ifs >> n;
			for(int i=0; i<n; i++) {
				metaInfo.clear();
				ifs >> elmID >> elmType >> nParam;
				for(int j=0; j<nParam; j++) {
					int a;
					ifs >> a;
					metaInfo.push_back(a);
				}
				if(elmType == 15) { // points defined by one index
					int a;
					ifs >> a;
					ptsOut << a << " " ;
					for(int j=0; j<nParam; j++)
						ptsOut << metaInfo[j] << " ";
					ptsOut << endl;
				} else if(elmType == 1) { // lines defined by two indices
					int a,b;
					ifs >> a >> b;
					lineOut << a << " " << b << " ";
					for(int j=0; j<nParam; j++)
						lineOut << metaInfo[j] << " ";
					lineOut << endl;
				} else if(elmType == 2) { // triangles defined by three corners
					int a,b,c;
					ifs >> a >> b >> c;
					triOut << a << " " << b << " " << c << " ";
					for(int j=0; j<nParam; j++)
						triOut << metaInfo[j] << " ";
					triOut << endl;
				} else if(elmType == 4) { // tetrahedrals defined by four corners
					int a,b,c,d;
					ifs >> a >> b >> c >> d;
					tetrOut << a << " " << b << " " << c << " " << d << " ";
					for(int j=0; j<nParam; j++)
						tetrOut << metaInfo[j] << " ";
					tetrOut << endl;
				} else {
					cerr << "Unsopported element code detected: " << elmType << ". Ignoring element and continuing\n";
					getline(ifs, inLine); // eat up the rest of the line and don't use it for anything
				}
			}
			getline(ifs, inLine); // should be just a endline character
			getline(ifs, inLine); // should be $EndElements, but don't bother checking
		} else if(inLine.compare("$PhysicalNames") == 0) {
			ifs >> n;
      for(int i=0; i<n; i++) {
        int dim, numb;
        ifs >> dim >> numb;
			  getline(ifs, inLine);
        cout << inLine << " (" << dim << "D) has ID tag: " << numb << endl;
      }
			getline(ifs, inLine); // should be $EndPhysicalNames, but don't bother checking
		} else if(inLine.compare("") == 0) {
			isDone = true;
		} else {
			cerr << "Unknown block detected. Writing all information gathered so far and terminating...\n";
			isDone = true;
		}
	}

	for(int i=0; i<nOutputFiles; i++)
		cout << "File " << *allFileName[i] << " written successfully" << endl;

	// close all files and exit cleanly
	nodeOut.close();
	ptsOut.close();
	lineOut.close();
	triOut.close();
	tetrOut.close();
	ifs.close();

}
