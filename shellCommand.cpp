#include "shellCommand.h"
#include <unistd.h>

IBD_SD::CShellCommand::CShellCommand(void)
{

}

IBD_SD::CShellCommand::~CShellCommand(void)
{

}

bool IBD_SD::CShellCommand::getformatList(std::string &cmd, std::vector<std::string> &fileList)
{
   //  std::cout<<"cmd"<<cmd<<std::endl;
    FILE *fpin = popen(cmd.c_str(), "r");//"../test/"

    if(NULL == fpin)
    {
	std::cout<<"fail to open file"<<std::endl;
        return false;
    }
    const size_t maxLine = 1000;
    char result[maxLine];
    while(1)
    {
        if( NULL == fgets(result, maxLine, fpin))
            break;
        std::string tmpString = result;
        if(tmpString[tmpString.size() - 1] == 10)
        {

            fileList.push_back(tmpString.substr(0, tmpString.size() - 1));
        }
        else{

            fileList.push_back(tmpString);
	}
    }
    if(0 != pclose(fpin))
    {
        return false;
    }
    return true;

}

bool IBD_SD::CShellCommand::getformatCount(std::string &cmd, int& count)
{
    FILE *fpin = popen(cmd.c_str(), "r");//"../test/"
    if(NULL == fpin)
    {
	std::cout<<"fail to open file"<<std::endl;
        return false;
    }
    const size_t maxLine = 1000;
    char result[maxLine];
    while(1)
    {
        if( NULL == fgets(result, maxLine, fpin))
            break;
        count++;
    }
    if(0 != pclose(fpin))
    {
        return false;
    }
    return true;
}


bool IBD_SD::CShellCommand::MakeDirection(std::string mkdir)
{
    std::string cmd="mkdir -p "+mkdir;
    FILE* fpin= popen(cmd.c_str(),"r");
    if(NULL==fpin)
       return false;
    if(0!=pclose(fpin))
       return false;
     return true;
}

void IBD_SD::CShellCommand::CreateProject(std::string project_path)
{
    std::string cmd1="mkdir -p "+project_path+"/ImageProcess/DeepImage/disparity_show/";
    std::string cmd2="mkdir -p "+project_path+"/ImageProcess/DetectImage/img_result/";
    std::string cmd3="mkdir -p "+project_path+"/ImageProcess/DetectImage/txt_result/";
    std::string cmd4="mkdir -p "+project_path+"/ImageProcess/3DResconstruct/";
    std::string cmd5="mkdir -p "+project_path+"/ImageProcess/3DResconstruct/Sign_Local3D/";

    std::string cmd6="mkdir -p "+project_path+"/ImageProcess/RectifyImage/left";
    std::string cmd7="mkdir -p "+project_path+"/ImageProcess/RectifyImage/right";

    std::string cmd8="mkdir -p "+project_path+"/ImageProcess/SemiImage/";
    std::string cmd9="mkdir -p "+project_path+"/ImageProcess/3DResconstruct/Lane_Local3D/";


    std::string cmd10="mkdir -p "+project_path+"/IMUProcess/";
    std::string cmd11="mkdir -p "+project_path+"/Rawdata/image/discardImg/";
    std::string cmd12="mkdir -p "+project_path+"/Rawdata/json_log/";

    std::string cmd13="mkdir -p "+project_path+"/VectorProcess/";
    std::string cmd14="mkdir -p "+project_path+"/ImageProcess/RectifyImage/Q";

    std::string cmd15="mkdir -p "+project_path+"/Log/";
    std::string cmd16="mkdir -p "+project_path+"/ImageProcess/3DResconstruct/Lane_Local3D/left/";

    std::string cmd17="mkdir -p "+project_path+"/ImageProcess/3DResconstruct/Lane_Local3D/right/";
    std::string cmd18="mkdir -p "+project_path+"/ImageProcess/SemiImage/left/";

    std::string cmd19="mkdir -p "+project_path+"/ImageProcess/SemiImage/right/";
    std::string cmd20="mkdir -p "+project_path+"/ImageProcess/SemiImage/lane_dashed/";

    std::string cmd21="mkdir -p "+project_path+"/ImageProcess/3DResconstruct/Lane_Local3D/left_deep/";
    std::string cmd22="mkdir -p "+project_path+"/VectorProcess/obj/";


    FILE* fpin1=popen(cmd1.c_str(),"r");
    pclose(fpin1);
    FILE* fpin2= popen(cmd2.c_str(),"r");
    pclose(fpin2);
    FILE* fpin3=popen(cmd3.c_str(),"r");
    pclose(fpin3);
    FILE* fpin4=popen(cmd4.c_str(),"r");
    pclose(fpin4);
    FILE* fpin5=popen(cmd5.c_str(),"r");
    pclose(fpin5);
    FILE* fpin6=popen(cmd6.c_str(),"r");
    pclose(fpin6);
    FILE* fpin7=popen(cmd7.c_str(),"r");
    pclose(fpin7);
    FILE* fpin8=popen(cmd8.c_str(),"r");
    pclose(fpin8);
    FILE* fpin9=popen(cmd9.c_str(),"r");
    pclose(fpin9);
    FILE* fpin10=popen(cmd10.c_str(),"r");
    pclose(fpin10);
    FILE* fpin11=popen(cmd11.c_str(),"r");
    pclose(fpin11);
    FILE* fpin12=popen(cmd12.c_str(),"r");
    pclose(fpin12);
    FILE* fpin13=popen(cmd13.c_str(),"r");
    pclose(fpin13);
    FILE* fpin14=popen(cmd14.c_str(),"r");
    pclose(fpin14);
    FILE* fpin15=popen(cmd15.c_str(),"r");
    pclose(fpin15);
    FILE* fpin16=popen(cmd16.c_str(),"r");
    pclose(fpin16);
    FILE* fpin17=popen(cmd17.c_str(),"r");
    pclose(fpin17);
    FILE* fpin18=popen(cmd18.c_str(),"r");
    pclose(fpin18);
    FILE* fpin19=popen(cmd19.c_str(),"r");
    pclose(fpin19);
    FILE* fpin20=popen(cmd20.c_str(),"r");
    pclose(fpin20);
    FILE* fpin21=popen(cmd21.c_str(),"r");
    pclose(fpin21);
    FILE* fpin22=popen(cmd22.c_str(),"r");
    pclose(fpin22);

    
}

//copyfile ---with absolute path
bool IBD_SD::CShellCommand::CopyFileToDstDir(std::string copyfile,std::string DstDirection)
{
    std::string cmd="cp "+copyfile+" "+DstDirection+"/";
    FILE* fpin=popen(cmd.c_str(),"r");
    if(NULL==fpin)
	fpin=popen(cmd.c_str(),"r");
    if(NULL==fpin)
	return false;
    if(!0==pclose(fpin))
	return false;
    return true;
}
//delet file
bool IBD_SD::CShellCommand::Deletefile(std::string deletefile)
{
    int flag=remove(deletefile.c_str());
    if(flag)
	return true;
    if(!flag)
	return false;
}
bool IBD_SD::CShellCommand::CutfileToDstdir(std::string srcFile,std::string dstDir)
{
    if(CopyFileToDstDir(srcFile,dstDir))
    {
	if(Deletefile(srcFile))
	    return true;
    }
    else
	return false;

}


bool IBD_SD::CShellCommand::CopyDirFile(std::string SrcDdir,std::string DstDir)
{
    std::string cmd="cp -r "+SrcDdir+"/. "+DstDir;
    LOG(INFO)<< "shellCommand_cpp : CopyDirFile : in...";

    FILE* fpin=popen(cmd.c_str(),"r");
    if(NULL==fpin)
       fpin=popen(cmd.c_str(),"r");
    if(NULL==fpin)
        return false;
    if(!0==pclose(fpin))
        return false;
    return true;

}

bool IBD_SD::CShellCommand::CopyRawdata(std::string Rawdata_path,std::string project_path)
{
    std::string src_file=Rawdata_path+"/image/";
    std::string dst_file=project_path+"/Rawdata/image/";
    if(!CopyDirFile(src_file,dst_file))
        return false;
    sleep(1);
    src_file=Rawdata_path+"/json_log/";
    dst_file=project_path+"/Rawdata/json_log/";
    if(!CopyDirFile(src_file,dst_file))
    {
        LOG(ERROR)<< "shellCommand_cpp : CopyRawdata :  failed !";
        return false;
    }
    else
    {
        LOG(INFO)<< "shellCommand_cpp : CopyRawdata :  successed !";
        return true;
    }
}











