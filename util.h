#ifndef DNECHO_UTIL_H
#define DNECHO_UTIL_H

typedef unsigned char BYTE;

//////////////////////////////////////////////////////////////////////////
//int KeyGrow(unsigned char * p, int w, int h, int* typeImg)
//函数用途：
//			查找图像连通区域，并记录连通区域相关信息，储存在Region region[REGION_NUM]中
//			返回找到的连通区域个数
//参数说明：
//p:	传入图像（opencv）的ImageData，要求该图像为二值图像，并且目标区域为黑点（像素值为0）
//w:	传入图像的宽
//h:	传入图像的高
//typeImg:
//int* typeImg=(int*)malloc(sizeof(int)*h*w);
//for(int i=0;i<h;i++)
//{
//	for(int j=0;j<w;j++)
//	{
//		typeImg[i*w+j]=-1;
//	}
//}
//////////////////////////////////////////////////////////////////////////


typedef struct Region
{
	int PointNum;
	int left;
	int right;
	int top;
	int bottom;
	int FusedN;
	int W;
	int H;

	int LeftTop;
	int RightTop;

	bool IsOK;
}Region;

#define MaxPointNum 20000
#define REGION_NUM 3000
extern Region region[REGION_NUM];

typedef struct CPoint
{
	int x;
	int y;
}CPoint;

typedef struct RectCenter{
	double x;
	double y;
}RectCenter;

int KeyGrow(unsigned char * p, int w, int h, int* typeImg);


//////////////////////////////////////////////////////////////////////////
//IplImage * RotateGrayImg(IplImage * src, double angle)
//函数用途：
//			将图片src旋转angle角度，顺时间方向转动为正角度，返回旋转后图像
//参数说明：
//src:		输入图像src，要求src是灰度图像
//angle:	需要旋转的角度，要求角度范围为[-90，90]
//////////////////////////////////////////////////////////////////////////
#define PI 3.1415926535
#define RADIAN(angle) ((angle)*PI/180.0) //角度到弧度转化的宏
IplImage * RotateGrayImg(IplImage * src, double angle);

//////////////////////////////////////////////////////////////////////////
//IplImage* RotateImage_CV(IplImage* img,float degree)
//函数用途：
//			将图片img旋转degree角度，顺时间方向转动为负角度，返回旋转后图像
//参数说明：
//src:		输入图像img
//angle:	需要旋转的角度
//////////////////////////////////////////////////////////////////////////
IplImage* RotateImage_CV(IplImage* img,float degree);

//////////////////////////////////////////////////////////////////////////
//void BiFilter(IplImage * src, IplImage * des,int halfL,float delta1,float delta2)
//函数用途：
//			双边滤波
//参数说明：
//src：		源图像
//des:		滤波后图像
//halfL:	双边滤波考虑的几何距离
//delta1:	距离高斯函数的标准差
//delta2:	像素高斯函数的标准差
//备注：	opencv的cvSmooth函数有双边滤波功能（CV_BILATERAL），双边滤波速度很慢，谨慎使用
//////////////////////////////////////////////////////////////////////////
void BiFilter(IplImage * src, IplImage * des,int halfL,float delta1,float delta2);

//////////////////////////////////////////////////////////////////////////
//IplImage* filterSigleChannel(IplImage* pSrc, char chn, int RTh, int GTh, int BTh)
//函数用途：
//			提取单通道信息，并返回二值化结果
//参数说明：
//pSrc：	源图像
//chn:		需要提取的单通道，可选参数为'r','g','b'
//Rth:		红色通道阈值
//Gth:		绿色通道阈值
//Bth:		蓝色通道阈值
//////////////////////////////////////////////////////////////////////////
IplImage* filterSigleChannel(IplImage* pSrc, char chn, int RTh, int GTh, int BTh);

//////////////////////////////////////////////////////////////////////////
//void myShowImg(char* win_name, IplImage* pSrc);
//函数用途：
//			快捷显示图片
//参数说明：
//win_name:	窗口名字
//pSrc：	源图像
//////////////////////////////////////////////////////////////////////////
void myShowImg(char* win_name, IplImage* pSrc);

//////////////////////////////////////////////////////////////////////////
//void MyGaussian(unsigned char *pSrcData, unsigned char *pDesData, int Width, int Height);
//函数用途：
//			高斯滤波
//参数说明：
//pSrcData:	源图像像素数据（内部已经根据Width处理了widthStep）
//pDesData：滤波后图像像素数据
//Width：	源图像宽度
//Height:	源图像高度
//////////////////////////////////////////////////////////////////////////
void MyGaussian(unsigned char *pSrcData, unsigned char *pDesData, int Width, int Height);

//////////////////////////////////////////////////////////////////////////
//void adaptiveThreshold_C(unsigned char* input, int IMAGE_WIDTH, int IMAGE_HEIGHT, int IMAGE_WIDESTEP, unsigned char* bin);
//函数用途：
//			自适应阈值分割
//参数说明：
//pSrcData:	源图像像素数据
//pDesData：滤波后图像像素数据
//IMAGE_WIDTH：	源图像宽度
//IMAGE_HEIGHT:	源图像高度
//IMAGE_WIDESTEP:	源图像步长
//S:	自适应阈值分割半径
//////////////////////////////////////////////////////////////////////////
void adaptiveThreshold_C(unsigned char* pSrcData, int IMAGE_WIDTH, int IMAGE_HEIGHT, int IMAGE_WIDESTEP, unsigned char* pDesData, int S=15, double T=0.15);

//////////////////////////////////////////////////////////////////////////
//void testImg(unsigned char* pData,int width, int height, int widthStep, char* win_name);
//函数用途：
//			显示pData对应的图像
//参数说明：
//pData:	图像像素数据
//width：	图像宽度
//height：	图像高度
//widthStep:	图像步长
//win_name:	显示窗口名字
//////////////////////////////////////////////////////////////////////////
void testImg(unsigned char* pData,int width, int height, int widthStep, char* win_name);

//////////////////////////////////////////////////////////////////////////
//函数用途：
//			获得指定目录下指定后缀的文件数量以及文件路径
//参数说明：
//dir:		文件目录最后不要求为\\，如"E:\\face"
//filename:	文件名字，包含目录，可直接用于打开文件
//filenum:	希望从目录中得到的文件数量
//suffix:	希望得到文件的后缀，如jpg，bmp（注意不用加.），默认为空，为空时不指定文件后缀名
//函数返回值：
//			表示实际得到的文件数量
//////////////////////////////////////////////////////////////////////////
int GetFileNameFromDir(char* _dir, char** filename, int filenum, char* suffix);

//////////////////////////////////////////////////////////////////////////
//函数用途：
//				将RGB图像转换为YCrCb图像
//参数说明：
//pSrcData：	RGB源图像，图像按opencv读入图像时的存贮格式，即BGR顺序
//width:		图像宽度
//height：		图像高度
//widthstep:	图像步长
//Y,Cr,Cb:		转换后Y,Cb,Cr通道数值
//////////////////////////////////////////////////////////////////////////
void myRGB2YCrCb(unsigned char* pSrcData,int width, int height, int widthstep, uchar* Y, uchar* Cr, uchar* Cb);


//////////////////////////////////////////////////////////////////////////
//函数用途：
//				用柱状图展示源图像的直方图
//参数说明：
//pSrc：		源图像（必须为单通道，深度depth为8）
//win_name:		显示图像的标题
//////////////////////////////////////////////////////////////////////////
void showHist(IplImage* pSrc, char* win_name);
#endif