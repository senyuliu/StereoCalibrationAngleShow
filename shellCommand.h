#ifndef __SHELLCOMMAND_H__
#define __SHELLCOMMAND_H__

#include<iostream>
#include<vector>
#include<stdio.h>
#include<stdlib.h>
#include<glog/logging.h>

using namespace std;

namespace IBD_SD
{

class CShellCommand
{
public: 

    CShellCommand(void);
    ~CShellCommand(void);

public: 

    bool getformatList(std::string &cmd, std::vector<std::string> &fileList);

    bool MakeDirection(std::string mkdir);
    bool getformatCount(std::string &cmd, int& count);

    void CreateProject(std::string project_path);
    bool CopyFileToDstDir(std::string copyfile,std::string DstDirection);

    bool Deletefile(std::string deletefile);
    bool CutfileToDstdir(std::string srcFile,std::string dstDir);

    bool CopyDirFile(std::string SrcDdir,std::string DstDir);
    bool CopyRawdata(std::string Rawdata_path,std::string project_path);

};

}

#endif 
