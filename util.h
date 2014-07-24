#ifndef DNECHO_UTIL_H
#define DNECHO_UTIL_H

#define bool int
#define false 0
#define true 1
typedef unsigned char BYTE;

//////////////////////////////////////////////////////////////////////////
//int KeyGrow(unsigned char * p, int w, int h)
//函数用途：
//			查找图像连通区域，并记录连通区域相关信息，储存在Region region[REGION_NUM]中
//			返回找到的连通区域个数
//参数说明：
//p:	传入图像（opencv）的ImageData，要求该图像为二值图像，并且目标区域为黑点（像素值为0）
//w:	传入图像的宽
//h:	传入图像的高
//////////////////////////////////////////////////////////////////////////
#define MaxPointNum 20000
#define REGION_NUM 1000
extern Region region[REGION_NUM];

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

typedef struct CPoint
{
	int x;
	int y;
}CPoint;

typedef struct RectCenter{
	double x;
	double y;
}RectCenter;

int KeyGrow(unsigned char * p, int w, int h);


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


#endif