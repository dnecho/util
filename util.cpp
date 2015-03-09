#include "stdafx.h"
#include "util.h"

//////////////////////////////////////////////////////////////////////////
//int KeyGrow(unsigned char * p, int w, int h)
//////////////////////////////////////////////////////////////////////////
Region region[REGION_NUM];
int KeyGrow(unsigned char * p, int w, int h, int* typeImg)
{
	int LineSize;
	unsigned char DealPixel;
	int * imgcopy;
	CPoint keytype[MaxPointNum];
	int N,nLEFT,nRight,nTOP,nBOTTOM;
	int i,j,n,typeNum;
	int x,y;
	
	if(p==0)
	{
		return 0;
	}

	//根据具体情况判断是否需要进行四字节对齐
	LineSize=(w+3)/4*4;
	DealPixel=0;

	imgcopy=(int*)malloc(sizeof(int)*LineSize*h);
	memset(imgcopy,0,sizeof(int)*LineSize*h);


	//类别的小标
	N=-1;
	//int nLEFT,nRight,nTOP,nBOTTOM;
	for (j=0;j<h;j++)
	{
		for (i=0;i<w;i++)
		{
			//数组下标
			n=-1;
			//当前白点的数目
			typeNum=0;
			//如果当前的是白点，并且还没有处理过，将其作为一个种子点，向四周扩散
			if(imgcopy[j*LineSize+i]==0 && p[j*LineSize+i]==DealPixel )
			{
				typeImg[j*w+i]=N+1;

				n++;
				keytype[n].x=i;
				keytype[n].y=j;
				typeNum++;

				imgcopy[j*LineSize+i]=1;		
			}
			else
				continue;

			nLEFT=LineSize;
			nRight=0;
			nTOP=h;
			nBOTTOM=0;
			while (n>=0 && n<MaxPointNum)
			{
				//数组末尾的值相当于POP出来
				x=keytype[n].x;
				y=keytype[n].y;

				nLEFT=nLEFT<x?nLEFT:x;
				nRight=nRight>x?nRight:x;
				nTOP=nTOP<y?nTOP:y;
				nBOTTOM=nBOTTOM>y?nBOTTOM:y;

				//下标减一
				n--;

				//zuo
				if(x-1>=0&&p[y*LineSize+x-1]==DealPixel&&imgcopy[y*LineSize+x-1]==0)
				{
					typeImg[y*w+x-1]=N+1;
					n++;

					if(n>=MaxPointNum)
					{
						break;
					}

					keytype[n].x=x-1;
					keytype[n].y=y;
					typeNum++;



					imgcopy[y*LineSize+x-1]=1;
				}
				//zuoshang
				if(x-1>=0&&y-1>=0&&p[(y-1)*LineSize+x-1]==DealPixel&&imgcopy[(y-1)*LineSize+x-1]==0)
				{
					typeImg[(y-1)*w+x-1]=N+1;
					n++;

					if(n>=MaxPointNum)
					{
						break;
					}

					keytype[n].x=x-1;
					keytype[n].y=y-1;
					typeNum++;



					imgcopy[(y-1)*LineSize+x-1]=1;
				}
				//youshang
				if(y-1>=0&&x+1<w&&p[(y-1)*LineSize+x+1]==DealPixel&&imgcopy[(y-1)*LineSize+x+1]==0)
				{
					typeImg[(y-1)*w+x+1]=N+1;
					n++;

					if(n>=MaxPointNum)
					{
						break;
					}

					keytype[n].x=x+1;
					keytype[n].y=y-1;
					typeNum++;



					imgcopy[(y-1)*LineSize+x+1]=1;
				}
				//you
				if(x+1<w&&p[y*LineSize+x+1]==DealPixel&&imgcopy[y*LineSize+x+1]==0)
				{
					typeImg[y*w+x+1]=N+1;
					n++;

					if(n>=MaxPointNum)
					{
						break;
					}

					keytype[n].x=x+1;
					keytype[n].y=y;
					typeNum++;



					imgcopy[y*LineSize+x+1]=1;
				}
				//shang
				if(y-1>=0&&p[(y-1)*LineSize+x]==DealPixel&&imgcopy[(y-1)*LineSize+x]==0)
				{
					typeImg[(y-1)*w+x]=N+1;
					n++;

					if(n>=MaxPointNum)
					{
						break;
					}

					keytype[n].x=x;
					keytype[n].y=y-1;
					typeNum++;



					imgcopy[(y-1)*LineSize+x]=1;
				}
				//zuoxia
				if(x-1>=0 && y+1<h&&p[(y+1)*LineSize+x-1]==DealPixel&&imgcopy[(y+1)*LineSize+x-1]==0)
				{
					typeImg[(y+1)*w+x-1]=N+1;
					n++;

					if(n>=MaxPointNum)
					{
						break;
					}

					keytype[n].x=x-1;
					keytype[n].y=y+1;
					typeNum++;



					imgcopy[(y+1)*LineSize+x-1]=1;
				}
				//xia
				if(y+1<h&&p[(y+1)*LineSize+x]==DealPixel&&imgcopy[(y+1)*LineSize+x]==0)
				{
					typeImg[(y+1)*w+x]=N+1;
					n++;


					if(n>=MaxPointNum)
					{
						break;
					}

					keytype[n].x=x;
					keytype[n].y=y+1;
					typeNum++;


					imgcopy[(y+1)*LineSize+x]=1;
				}
				//youxia
				if(y+1<h&&x+1<w&&p[(y+1)*LineSize+x+1]==DealPixel&&imgcopy[(y+1)*LineSize+x+1]==0)
				{
					typeImg[(y+1)*w+x+1]=N+1;
					n++;

					if(n>=MaxPointNum)
					{
						break;
					}

					keytype[n].x=x+1;
					keytype[n].y=y+1;
					typeNum++;



					imgcopy[(y+1)*LineSize+x+1]=1;
				}
			}
			//过滤后保存在数组当中的裂纹为预选的裂纹,后续再进行合并
			if(typeNum>0)
			{
				if (N>=REGION_NUM)
				{
					free(imgcopy);
					return REGION_NUM;
				}

				N++;
				region[N].PointNum=typeNum;
				region[N].left=nLEFT;
				region[N].right=nRight;
				region[N].top=nTOP;
				region[N].bottom=nBOTTOM;
				region[N].IsOK=true;
				region[N].FusedN=0;
				region[N].W=nRight-nLEFT+1;
				region[N].H=nBOTTOM-nTOP+1;

				region[N].LeftTop=nTOP;
				region[N].RightTop=nTOP;
			}
		}
	}

	free(imgcopy);

	return ++N;
}

//////////////////////////////////////////////////////////////////////////
//IplImage * RotateGrayImg(IplImage * src, double angle)
//////////////////////////////////////////////////////////////////////////
IplImage * RotateGrayImg(IplImage * src, double angle)
{
//	int m_Angle=angle;
	int width=src->width;
	int height=src->height;

	BYTE * lpTempPtr=NULL;
	BYTE * p_Copy=(BYTE *)src->imageData;

	double       SrcX1,SrcY1,SrcX2,SrcY2;
	double       SrcX3,SrcY3,SrcX4,SrcY4;
	double       DstX1,DstY1,DstX2,DstY2;
	double       DstX3,DstY3,DstX4,DstY4;
	double       num1,num2;

	double RotateAngle=(double)RADIAN(angle);
	//double RotateAngle=angle;

	double cosa=(double)cos((double)RotateAngle);
	double sina=(double)sin((double)RotateAngle);
	//原图的宽度和高度
	//   BYTE * lpSrc=p_Copy;
	int Wold=width;
	int Hold=height;
	//原图的四个角的坐标
	SrcX1=(double)(-0.5*Wold);
	SrcY1=(double)(0.5*Hold);
	SrcX2=(double)(0.5*Wold);
	SrcY2=(double)(0.5*Hold);
	SrcX3=(double)(-0.5*Wold);
	SrcY3=(double)(-0.5*Hold);
	SrcX4=(double)(0.5*Wold);
	SrcY4=(double)(-0.5*Hold);
	//新图四个角的坐标
	DstX1=cosa*SrcX1+sina*SrcY1;
	DstY1=-sina*SrcX1+cosa*SrcY1;
	DstX2=cosa*SrcX2+sina*SrcY2;
	DstY2=-sina*SrcX2+cosa*SrcY2;
	DstX3=cosa*SrcX3+sina*SrcY3;
	DstY3=-sina*SrcX3+cosa*SrcY3;
	DstX4=cosa*SrcX4+sina*SrcY4;
	DstY4=-sina*SrcX4+cosa*SrcY4;
	//计算新图的宽度，高度
	int Wnew = (int)(max(fabs(DstX4-DstX1), fabs(DstX3-DstX2))+0.5);
	int Hnew = (int)(max(fabs(DstY4-DstY1), fabs(DstY3-DstY2))+0.5);
	//计算矩阵(2.9)中的两个常数，这样不用以后每次都计算了
	num1=(double)( -0.5*Wnew*cosa-0.5*Hnew*sina+0.5*Wold);
	num2=(double)(0.5*Wnew*sina-0.5*Hnew*cosa+0.5*Hold);

	//将新的缓冲区中的每个字节都填成255，这样以后未处理的象素就是白色
	int DstBufSize=Wnew*Hnew;

	/*if(lpTempPtr==NULL)
	{
		lpTempPtr=new BYTE[Wnew*Hnew];
		memset(lpTempPtr,(BYTE)255,Wnew*Hnew);
	}
	else
	{
		delete []lpTempPtr;
		lpTempPtr=NULL;

		lpTempPtr=new BYTE[Wnew*Hnew];
		memset(lpTempPtr,(BYTE)255,Wnew*Hnew);
	}*/

	IplImage *Des=cvCreateImage(cvSize(Wnew,Hnew),src->depth,src->nChannels);
	lpTempPtr=(BYTE *)Des->imageData;
	memset(lpTempPtr,255,Des->widthStep*Des->height);

	int        x0,y0,x1,y1;
	int nChannels = src->nChannels;
	for(y1=0;y1<Hnew;y1++)
	{
		for(x1=0;x1<Wnew;x1++)
		{
			//x0,y0为对应的原图上的坐标
			//for(int iChannel = 0; iChannel < nChannels; iChannel++)
			//{
			//	x0= (int)((x1*nChannels+iChannel)*cosa+y1*sina+num1);
			//	y0= (int)(-1.0f*(x1*nChannels+iChannel)*sina+y1*cosa+num2);
			//	if( (x0>=0) && (x0<src->widthStep) && (y0>=0) && (y0<Hold))   //在原图范围内
			//	{
			//		lpTempPtr[y1*Des->widthStep+x1*nChannels+iChannel]=p_Copy[y0*src->widthStep+x0];
			//	}
			//}

			x0= (int)(x1*cosa+y1*sina+num1);
			y0= (int)(-1.0f*x1*sina+y1*cosa+num2);
			if( (x0>=0) && (x0<Wold) && (y0>=0) && (y0<Hold))   //在原图范围内
			{
				lpTempPtr[y1*Des->widthStep+x1]=p_Copy[y0*src->widthStep+x0];
			}
		}
	}

//	MyNamedWindow("rr",1);
//	MyShowImage("rr",Des);

	return Des;
}

//////////////////////////////////////////////////////////////////////////
//IplImage* RotateImage_CV(IplImage* img,float degree)
//////////////////////////////////////////////////////////////////////////
IplImage* RotateImage_CV(IplImage* img,float degree)
{  
	double angle = degree  * CV_PI / 180.; // 弧度    
	double a = sin(angle), b = cos(angle);   
	int width = img->width;    
	int height = img->height;    
	int width_rotate= int(height * fabs(a) + width * fabs(b));    
	int height_rotate=int(width * fabs(a) + height * fabs(b));    
	//旋转数组map  
	// [ m0  m1  m2 ] ===>  [ A11  A12   b1 ]  
	// [ m3  m4  m5 ] ===>  [ A21  A22   b2 ]  
	float map[6]={0};  
	CvMat map_matrix = cvMat(2, 3, CV_32F, map);    
	// 旋转中心  
	CvPoint2D32f center = cvPoint2D32f(width / 2, height / 2);    
	cv2DRotationMatrix(center, degree, 1.0, &map_matrix);    
	map[2] += (width_rotate - width) / 2;    
	map[5] += (height_rotate - height) / 2;    
	IplImage* img_rotate = cvCreateImage(cvSize(width_rotate, height_rotate), 8, img->nChannels);   
	//对图像做仿射变换  
	//CV_WARP_FILL_OUTLIERS - 填充所有输出图像的象素。  
	//如果部分象素落在输入图像的边界外，那么它们的值设定为 fillval.  
	//CV_WARP_INVERSE_MAP - 指定 map_matrix 是输出图像到输入图像的反变换，  
	cvWarpAffine( img,img_rotate, &map_matrix, CV_INTER_LINEAR | CV_WARP_FILL_OUTLIERS, cvScalarAll(0));    
	return img_rotate;  
}  

//////////////////////////////////////////////////////////////////////////
//void BiFilter(IplImage * src, IplImage * des,int halfL,float delta1,float delta2)
//////////////////////////////////////////////////////////////////////////
void BiFilter(IplImage * src, IplImage * des,int halfL,float delta1,float delta2)
{
	int i,j,m,n;
	int height=src->height;
	int width=src->width;
	BYTE * Data=(BYTE *)src->imageData;

	//高斯模糊的权重
	int L=2*halfL+1;
	float * gaussW=new float[L*L];
	for (i=-halfL;i<=halfL;i++)
	{
		for (j=-halfL;j<=halfL;j++)
		{
			gaussW[(i+halfL)*L+(j+halfL)]=(float)(1/sqrt(2*3.1415926*delta1)*pow(2.718281828,-1.0*(i*i+j*j)/(2*delta1)));
		}
	}

	float * gaussPix=new float[256];
	for (j=0;j<256;j++)
	{
		gaussPix[j]=(float)(1/sqrt(2*3.1415926*delta2)*pow(2.718281828,-1.0*j*j/(2*delta2)));
	}

	//开始模糊
	float pixV,weightsum;
	for (i=0;i<height;i++)
	{
		for (j=0;j<width;j++)
		{
			for (int mm=0;mm<src->nChannels;mm++)
			{
				weightsum=0;
				pixV=0;
				for (m=-halfL;m<=halfL;m++)
				{
					for (n=-halfL;n<=halfL;n++)
					{
						if (i+m<0 || i+m>height-1 || j+n<0 ||j+n>width-1)
						{

						}
						else
						{
							int pos=abs(Data[(i+m)*src->widthStep+(j+n)*src->nChannels+mm]-Data[i*src->widthStep+j*src->nChannels+mm]);
							float W= (gaussW[(m+halfL)*L+(n+halfL)]*gaussPix[pos]);
							weightsum+=W;
							pixV=pixV+W*Data[(i+m)*src->widthStep+(j+n)*src->nChannels+mm];
						}
					}
				}

				des->imageData[i*des->widthStep+j*des->nChannels+mm]=(char)((pixV+1e-6)/(weightsum+1e-6));
			}
		}
	}

	delete []gaussPix;
	delete []gaussW;
}

//////////////////////////////////////////////////////////////////////////
//IplImage* filterSigleChannel(IplImage* pSrc, char chn, int RTh, int GTh, int BTh)
//////////////////////////////////////////////////////////////////////////
IplImage* filterSigleChannel(IplImage* pSrc, char chn, int RTh, int GTh, int BTh)
{
	IplImage *pSrcR = cvCreateImage(cvGetSize(pSrc), pSrc->depth, 1);
	IplImage *pSrcG = cvCreateImage(cvGetSize(pSrc), pSrc->depth, 1);
	IplImage *pSrcB = cvCreateImage(cvGetSize(pSrc), pSrc->depth, 1);
	cvSplit(pSrc, pSrcB, pSrcG, pSrcR, NULL);

	unsigned char* pData = (unsigned char*)pSrc->imageData;
	unsigned char rgb[3];
	int th[3];
	IplImage *pTemp = cvCreateImage(cvGetSize(pSrc), pSrc->depth, 1);
	unsigned char* pTempData = (unsigned char*)pTemp->imageData;

	for(int i = 0; i < pSrc->height; i++)
	{
		for(int j = 0; j < pSrc->width; j++)
		{
			int pos = i*pSrc->widthStep+j*3;
			switch(chn){
				case 'r':
					th[2]=RTh; th[1]=GTh; th[0]=BTh;
					rgb[2] = pData[pos+2];//R
					rgb[1] = pData[pos+1];//G
					rgb[0] = pData[pos];//B
					break;
				case 'g':
					th[2]=GTh; th[1]=BTh; th[0]=RTh;
					rgb[2] = pData[pos+1];//G
					rgb[1] = pData[pos];//B
					rgb[0] = pData[pos+2];//R
					break;
				case 'b':
					th[2]=BTh; th[1]=GTh; th[0]=RTh;
					rgb[2] = pData[pos];//B
					rgb[1] = pData[pos+1];//G
					rgb[0] = pData[pos+2];//R
					break;
				default:
					cvReleaseImage(&pSrcR);
					cvReleaseImage(&pSrcG);
					cvReleaseImage(&pSrcB);
					return pTemp;

			}

			//if( rgb[1]<gbT && rgb[0]<gbT && rgb[2]>rT)
			if( rgb[0]<th[0] && rgb[1]<th[1] && rgb[2]>th[2])
			{
				pTempData[i*(pTemp->widthStep)+j] = 0;
			}else
			{
				pTempData[i*(pTemp->widthStep)+j] = 255;
			}
		}
	}
	cvReleaseImage(&pSrcR);
	cvReleaseImage(&pSrcG);
	cvReleaseImage(&pSrcB);
	return pTemp;
}

//////////////////////////////////////////////////////////////////////////
//void myShowImg(char* win_name, IplImage* pSrc)
//////////////////////////////////////////////////////////////////////////
void myShowImg(char* win_name, IplImage* pSrc)
{
	cvNamedWindow(win_name, 0);
	cvShowImage(win_name, pSrc);
}


//////////////////////////////////////////////////////////////////////////
//void MyGaussian(unsigned char *pSrcData, unsigned char *pDesData, int Width, int Height)
//多通道图像高斯平滑
//////////////////////////////////////////////////////////////////////////
void MyGaussian(unsigned char *pSrcData, unsigned char *pDesData, int Width, int Height)
{
	int small_gaussian[4][4] =
	{
		{1024}, 
		{512, 256}, 
		{384 , 256, 64},
		{288, 224, 112, 32}
	};

	int i,j,k,m;
	//int Width=pSrc->width;
	//int Height=pSrc->height;
	int nCh=1;
	int widestep=(Width+3)/4*4;
	memcpy(pDesData, pSrcData, sizeof(unsigned char)*Height*widestep);
	unsigned char * Data=pSrcData;	
	unsigned char * sData;
	int SmoothRadio=3;
	int pos,rpos;
	int Temp;
	int *  pGauTepData=NULL;
	unsigned char * rfilter;

	//根据SMOOTHRADIO判断是哪个数组
	if (SmoothRadio==1)
	{
		pGauTepData=small_gaussian[1];
	}
	else if(SmoothRadio==2)
	{
		pGauTepData=small_gaussian[2];
	}
	else if(SmoothRadio==3)
	{
		pGauTepData=small_gaussian[3];
	}
	else
	{

	}

	//rfilter=(unsigned char*)malloc(sizeof(unsigned char)*Width*Height);
	rfilter = new unsigned char[sizeof(unsigned char)*widestep*Height];
	sData=rfilter;
	//行滤波
	pos=0;
	rpos=0;
	for (i=0;i<Height;i++)
	{
		for (j=3;j<Width-3;j++)
		{
			rpos=pos+j*nCh;
			for (m=0;m<nCh;m++)
			{
				Temp=Data[rpos+m]*pGauTepData[0];
				for (k=1;k<=SmoothRadio;k++)
				{
					Temp+=(Data[rpos+m-k*nCh]+Data[rpos+m+k*nCh])*pGauTepData[k];
				}

				sData[rpos+m]=Temp>>10;
			}
		}
		pos+=widestep;
	}

	sData=pDesData;
	Data=rfilter;
	pos=3*widestep;
	rpos=0;
	for (i=3;i<Height-3;i++)
	{
		for (j=3;j<Width-3;j++)
		{
			rpos=pos+j*nCh;
			for (m=0;m<nCh;m++)
			{
				Temp=Data[rpos+m]*pGauTepData[0];
				for (k=1;k<=SmoothRadio;k++)
				{
					Temp+=(Data[rpos+m-k*widestep]+Data[rpos+m+k*widestep])*pGauTepData[k];
				}

				sData[rpos+m]=Temp>>10;
			}
		}
		pos+=widestep;
	}

	delete []rfilter;
	//free(rfilter);
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
void adaptiveThreshold_C(unsigned char* pSrcData, int IMAGE_WIDTH, int IMAGE_HEIGHT, int IMAGE_WIDESTEP, unsigned char* pDesData, int S, double T)
{
	unsigned long* integralImg = 0;
	int i, j;
	long sum=0;
	int count=0;
	int index;
	int x1, y1, x2, y2;
	int s2 = S/2;
	//t_start = clock();
	// create the integral image
	integralImg = (unsigned long*)malloc(IMAGE_WIDESTEP*IMAGE_HEIGHT*sizeof(unsigned long));
	memset(integralImg,0,sizeof(unsigned long*)*IMAGE_WIDESTEP*IMAGE_HEIGHT);

	for (i=0; i<IMAGE_WIDTH; i++)
	{
		// reset this column sum
		sum = 0;

		for (j=0; j<IMAGE_HEIGHT; j++)
		{
			index = j*IMAGE_WIDESTEP+i;

			sum += pSrcData[index];
			if (i==0)
				integralImg[index] = sum;
			else
				integralImg[index] = integralImg[index-1] + sum;
		}
	}
	/*t_end = clock();
	printf("integralImg:%d\n", t_end-t_start);
	t_start = clock();*/
	// perform thresholding
	for (i=0; i<IMAGE_WIDTH; i++)
	{
		for (j=0; j<IMAGE_HEIGHT; j++)
		{
			index = j*IMAGE_WIDESTEP+i;

			// set the SxS region
			x1=i-s2; x2=i+s2;
			y1=j-s2; y2=j+s2;

			// check the border
			if (x1 < 0) x1 = 0;
			if (x2 >= IMAGE_WIDTH) x2 = IMAGE_WIDTH-1;
			if (y1 < 0) y1 = 0;
			if (y2 >= IMAGE_HEIGHT) y2 = IMAGE_HEIGHT-1;

			count = (x2-x1)*(y2-y1);

			// I(x,y)=s(x2,y2)-s(x1,y2)-s(x2,y1)+s(x1,x1)
			sum = integralImg[y2*IMAGE_WIDESTEP+x2] -
				integralImg[y1*IMAGE_WIDESTEP+x2] -
				integralImg[y2*IMAGE_WIDESTEP+x1] +
				integralImg[y1*IMAGE_WIDESTEP+x1];

			if ((long)(pSrcData[index]*count) < (long)(sum*(1.0-T)))
				pDesData[index] = 0;
			else
				pDesData[index] = 255;
		}
	}
	
	free (integralImg);
	//t_end = clock();
	//printf("Threshold:%d\n", t_end-t_start);
}


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
void testImg(unsigned char* pData,int width, int height, int widthStep, char* win_name)
{
	IplImage* pImg = cvCreateImageHeader(cvSize(width,height), 8, 1);
	cvSetData(pImg, pData, widthStep);
	cvNamedWindow(win_name, 0);
	cvShowImage(win_name, pImg);
	cvWaitKey();
	cvReleaseImageHeader(&pImg);
}

//////////////////////////////////////////////////////////////////////////
//int GetFileNameFromDir(char* _dir, char** filename, int filenum, char* suffix)
//////////////////////////////////////////////////////////////////////////
int GetFileNameFromDir(char* _dir, char** filename, int filenum, char* suffix)
{
	int num = 0;
	int len=strlen(_dir);
	char dir[255];
	if(dir[len - 1] != '\\')
		sprintf_s(dir, "%s\\",_dir);
	else
		strcpy_s(dir, _dir);
	_finddata_t file;
	long longf;
	char tep_char[255];

	//得到需要查找的路径及需要查找文件的类型（是否指定后缀名）
	if(suffix == NULL)//判断是否需要制定文件后缀
	{
		sprintf_s(tep_char, "%s*.*", dir);
	}
	else
	{		
		sprintf_s(tep_char, "%s*.%s", dir, suffix);
	}

	//_findfirst返回的是long型; long __cdecl _findfirst(const char *, struct _finddata_t *)
	if((longf = _findfirst(tep_char, &file))==-1l)
	{
		printf("文件没有找到!\n");
	}
	else
	{
		if (file.name[0] != '.')//当后缀是"null"时findfirst找到的是'.'；如果后缀有值是找到的是第一个对应后缀的文件，必须保存起来
		{
			sprintf_s(tep_char, "%s%s", dir, file.name);//将路径与文件名结合
			strcpy_s(filename[num], 255, tep_char);
			num++;
		}

		//int __cdecl _findnext(long, struct _finddata_t *);如果找到下个文件的名字成功的话就返回0,否则返回-1
		while( _findnext( longf, &file ) == 0 && num <filenum)
		{
			if (file.name[0] == '.' && file.name[1] == '.')
			{
				continue;
			}
			sprintf_s(tep_char, "%s%s", dir, file.name);//将路径与文件名结合
			strcpy_s(filename[num], 255, tep_char);
			num++;
		}
	}

	_findclose(longf);
	return num;
}

//////////////////////////////////////////////////////////////////////////
//void myRGB2YCrCb(unsigned char* pSrcData,int width, int height, int widthstep, uchar* Y, uchar* Cr, uchar* Cb)
//////////////////////////////////////////////////////////////////////////
void myRGB2YCrCb(unsigned char* pSrcData,int width, int height, int widthstep, uchar* Y, uchar* Cr, uchar* Cb)
{
	int R,G,B,tmp;
	unsigned char* Y_tmp = Y;
	unsigned char* Cb_tmp = Cb;
	unsigned char* Cr_tmp = Cr;
	for(int i=0; i<height; i++)
	{
		unsigned char* pRowData = pSrcData+i*widthstep;

		for(int j=0; j<width; j++)
		{
			B = *pRowData++;
			G = *pRowData++;
			R = *pRowData++;

			*Y_tmp++ = 0.299*R+0.587*G+0.114*B;
			*Cb_tmp++ = (-0.168*R - 0.3313*G + 0.5*B)+128;
			*Cr_tmp++ = (0.5*R -0.4187*G -0.0813*B)+128;
		}
	}
}