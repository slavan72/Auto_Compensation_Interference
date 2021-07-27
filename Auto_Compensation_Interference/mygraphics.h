# if !defined(__MYGRAPHICS)
#define __MYGRAPHICS

#include "windows.h"
#include "windowsx.h"


#define  LEFT_TEXT		0
#define  CENTER_TEXT	1
#define  RIGHT_TEXT		2
#define  BOTTOM_TEXT	0
#define  TOP_TEXT		2
#define  DEFAULT_FONT	0
#define  HORIZ_DIR		0
#define  VERT_DIR		1
#define  DETECT			0
#define  grOk			0
#define	 MAXCOLORS	   128


struct palettetype{
	unsigned char size;
	signed char colors[MAXCOLORS+1];
};

//conio.h
void window(int left, int top, int right, int bottom){}

void textattr(int attr){
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), attr);
}

void clrscr(){
	system("cls");
}

void gotoxy(int x, int y){
	COORD position = {x, y};
	SetConsoleCursorPosition(GetStdHandle(STD_OUTPUT_HANDLE), position);
}

//graphics.h
void settextstyle(int font, int direction, int charsize){}
void settextjustify(int horiz, int vert){}
void setpalette(int colornum, int color){}
void setcolor(int color){}
void setbkcolor(int color){}

void rectangle(int left, int top, int right, int bottom){
	HWND handle = GetConsoleWindow(); //HWND handle = FindWindow("ConsoleWindowClass", NULL);
	HDC hdc = GetDC(handle);

	Rectangle(hdc, left, top, right, bottom);
}
void bar(int left, int top, int right, int bottom){

	HWND handle = GetConsoleWindow(); //HWND handle = FindWindow("ConsoleWindowClass", NULL);
	HDC hdc = GetDC(handle);
}
void circle(int x, int y, int radius){}
void pieslice(int x, int y, int start, int end, int radius){}
void putpixel(int x, int y, int color){

	HWND handle = GetConsoleWindow();
	HDC hdc = GetDC(handle);

	SetPixel(hdc, x, y, color);
}

void outtextxy(int x, int y, char *str){

	HWND handle = GetConsoleWindow();
	HDC hdc = GetDC(handle);

	TextOutA(hdc, x, y, str, 2);
}
void setusercharsize(int multx, int divx, int multy, int divy){}
void setfillstyle(int patern, int color){}
void linerel(int deltax, int deltay){}
void getpalette(palettetype *palette){}
void moveto(int x, int y){}
void initgraph(int *graphdriver, int *graphmode, char *pathtodriver){}
void getaspectratio(int *xasp, int *yasp){}
void closegraph(){}
void cleardevice(){}

int textwidth(char *str){return sizeof(str);}
int textheight(char *str){return sizeof(str);}
int getmaxcolor(){
	int maxColor = 256;
	return maxColor;}
int getmaxx(){
	int maxX = 256;
	return maxX;}
int getmaxy(){
	int maxY;
	return maxY = 256;}
int graphresult(){
	int errorCode  = 0;
	return errorCode;
}
char* grapherrormsg(int errorCode){
	char* errCode = "ERROR";
	return errCode;
}




#endif __MYGRAPHICS