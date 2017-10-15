#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "ap.h"
#include "resource.h"
#include "resrc1.h"
#define dimmas 301   				//<-size of the picture + 1

 
LRESULT CALLBACK WindowFunc(HWND,UINT,WPARAM,LPARAM);
LRESULT CALLBACK WindowSon(HWND,UINT,WPARAM,LPARAM);
BOOL CALLBACK DialogProc(HWND hDlg, UINT uMsg, WPARAM wParam, LPARAM lParam);

double polynomial(int n, double ro, double theta);
long double findaver(void);
long double response(int kontakt, double r, double theta);
void ludecomposition(ap::real_2d_array& a,
     int m,
     int n,
     ap::integer_1d_array& pivots);
void ludecompositionunpacked(ap::real_2d_array a,
     int m,
     int n,
     ap::real_2d_array& l,
     ap::real_2d_array& u,
     ap::integer_1d_array& pivots);
bool solvesystemlu(const ap::real_2d_array& a,
     const ap::integer_1d_array& pivots,
     ap::real_1d_array b,
     int n,
     ap::real_1d_array& x);
bool solvesystem(ap::real_2d_array a,
	 ap::real_1d_array b,
	 int n,
	 ap::real_1d_array& x);

long double findnorm(void); 			// ищет норму для matrix для вывода на экран.

const int number_contacts=17;
const int step=10; 				// шаг выборки по каждой коорд. для МНК.

//////////// Глобальные переменные: начало //////////////
  char szWinName[] = "Мое Окно",
       WinSonName[]= "Окно для отрисовки";
  HWND hwnd,hwndson,hwndDlg=NULL;
  HMENU hMenu;
  HDC memdc,memdc1;
  PAINTSTRUCT paintstruct;
  HINSTANCE hInstance;
  LPSTR CommandLine, lpbuf="CheCK";
  unsigned long double PI;
  float tmp1, tmp2,tmp3;
  double nachs;
  double zernike[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double response[number_contacts][13]=
	{{0.00716,0.00008,-0.00007,-0.00023,-0.00001,0.00022,-0.00008,0.00006,-0.00011,-0.00003,0.00018,-0.00003,-0.00004},
	{0.00106,0.000755,0.000055,-0.00001,0.00079,0.000045,-0.00001,0.00014,-0.000295,-0.000025,0.000035,-0.000105,-0.000055},
	{0.00102,-0.00009,-0.000795,0.00066,0.000575,0.000065,0.000055,-0.00004,-0.00001,0.00028,-0.000105,-0.00006,-0.000055},
	{0.00095,-0.000845,-0.0002,0.000965,0.000065,0.00009,-0.000095,-0.00001,0.00032,0.000055,-0.00015,0.000025,-0.000065},
	{0.000995,-0.000125,0.00071,0.000735,-0.00056,0.000045,0.00004,0.00011,0.0001,-0.00025,-0.00009,0.000115,-0.000055},
	{0.00108,0.000705,0.00025,0.00018,-0.00095,-0.000005,0.00006,-0.000005,-0.0002,-0.00008,-0.000005,0.0002,-0.00005},
	{0.00105,0.000165,-0.00063,-0.000595,-0.00077,-0.000005,-0.00011,0.00006,-0.000045,0.000325,0.00011,0.00018,-0.000055},
	{0.000985,-0.000615,-0.00017,-0.001,-0.000185,0.000015,0.00013,0.00011,0.00033,0.00009,0.00023,0.000015,-0.00008},
	{0.00108,-0.00008,0.000635,-0.000815,0.00054,-0.000035,0.000005,-0.00005,0.00005,-0.00028,0.000145,-0.000135,-0.00007},
	{0.000325,0.001175,-0.00086,0.00009,0.00021,0.0001,0.00028,0.00025,-0.000265,0.000215,0.000055,0.00013,0.00001},
	{0.00031,-0.00091,-0.001095,0.00031,-0.00004,0.00012,-0.000115,-0.000385,0.00012,0.00027,0.000105,0.00013,0.000025},
	{0.000355,-0.000995,0.00073,0.00028,-0.00013,0.000135,-0.00022,0.00039,0.000245,-0.000125,0.000125,-0.00001,-0.00002},
	{0.000335,0.000875,0.001025,0.00011,-0.000225,0.000115,0.00046,-0.0002,-0.00023,-0.000275,0.000055,-0.000145,0.000005},
	{0.00029,0.001095,-0.000855,-0.000055,-0.000245,0.00009,-0.00039,-0.00024,-0.000285,0.00018,-0.000045,-0.000115,0.00002},
	{0.000305,-0.00079,-0.00107,-0.00021,-0.000265,0.000115,0.00014,0.000445,0.000205,0.000235,-0.000105,-0.00004,-0.000025},
	{0.00027,-0.00102,0.000865,-0.0002,-0.000025,0.000115,0.000255,-0.000425,0.00016,-0.000165,-0.000115,0.00008,0.00001},
	{0.00028,0.000885,0.000995,-0.000105,0.000265,0.0001,-0.0005,0.000085,-0.00027,-0.000175,-0.00005,0.000135,0.00002}};
  
  double norma=1;
  long global_t=10000;	   			//value for my timer
  BOOL eraser=FALSE,eraser1=FALSE;
  int ns;				   	//для решения системы:
  ap::real_2d_array square_matrix; 		//матрица коэффициентов [1..N][1..N]
  ap::real_1d_array column; 			//столбец свободных членов
  ap::real_1d_array voltages; 			//искомый столбец напряжений
  long double matrix[dimmas][dimmas],matrix1[dimmas][dimmas];
  char mybuf[1000]="proba";

  //////////// Глобальные переменные: конец ///////////////

  int WINAPI WinMain (HINSTANCE hThisInst, 
  			HINSTANCE hPrevInst,
  			LPSTR lpszArgs,
  			int nWinMode)
  {
  	MSG msg;
  	WNDCLASS wcl,wclson;

	hInstance=GetModuleHandle(NULL);
    
  // Определить класс окна 
 
  	wcl.hInstance=hThisInst;
  	wcl.lpszClassName=szWinName; 			//имя класса окна
  	wcl.lpfnWndProc=WindowFunc; 			//функция окна
  	wcl.style=0;                			// стиль по умолчанию
  	wcl.hIcon=LoadIcon(NULL,IDI_EXCLAMATION); 	//стандартная иконка
    	wcl.hCursor=LoadCursor(NULL,IDC_ARROW);
   	wcl.lpszMenuName= (LPSTR)MYMENU; 		//без меню  
    	wcl.cbClsExtra=0;
    	wcl.cbWndExtra=0;
    	// Заполнить окно цветом
    	wcl.hbrBackground=(HBRUSH) GetStockObject(BLACK_BRUSH);
    
	wclson=wcl;
	wclson.lpszClassName=WinSonName;
    	wclson.lpfnWndProc=WindowSon;
  	wclson.lpszMenuName=NULL;

    if ((!RegisterClass (&wcl))||(!RegisterClass (&wclson))) return 0;
    
    // Создать окно
    hwnd = CreateWindow(szWinName,"Курсовая",
    	   WS_OVERLAPPEDWINDOW,
	  	   100, 		// горизонтальное положение окна
	  	   100, 		// вертикальное положение окна	
	  	   dimmas+6, 		// ширина окна	
	  	   dimmas+25, 		// высота окна	 
	  	   HWND_DESKTOP,
	  	   NULL,
	  	   hThisInst,
	  	   NULL);

	hwndson = CreateWindow(WinSonName, "Окно для отрисовки",
    	   WS_CAPTION,
	  	   450, 		// горизонтальное положение окна
	  	   100, 		// вертикальное положение окна	
	  	   dimmas+6, 		// ширина окна	
	  	   dimmas+25, 		// высота окна	 
	  	   HWND_DESKTOP,
	  	   NULL,
	  	   hThisInst,
	  	   NULL);

	// Показать окно и нарисовать содержимое
	ShowWindow(hwndson,nWinMode);
    	ShowWindow(hwnd,nWinMode);
	UpdateWindow (hwndson);
	UpdateWindow (hwnd);
					
	// Цикл обработки сообщений
	while (GetMessage (&msg, NULL, 0, 0))
	{
		if ((hwndDlg!=0)&&(IsDialogMessage(hwndDlg,&msg))) continue;
		TranslateMessage (&msg);
		DispatchMessage (&msg);
	}

	DeleteDC (memdc);
	return msg.wParam;
  }	                  
  

LRESULT CALLBACK WindowFunc(HWND hwnd,UINT message,WPARAM wParam,LPARAM lParam)
  {
	int i,j,x,y,k,s;
	double theta,ro,theta,r;
	HDC hdc;
	HBITMAP hbit;
	HBRUSH hbrush;
	static int maxX, maxY;  

  	switch (message)
  	{
  	
		case WM_CREATE:
		PI=acos(-1);
		/*     	GetPrivateProfileString("main","global_t","error", &mybuf,20,iniad);
			global_t = atof(&mybuf);
			GetPrivateProfileString("main","global_a","error", &mybuf,20,iniad);
			global_a = atof(&mybuf);
		*/	eraser=TRUE;
			eraser1=TRUE;
		    	maxX = GetSystemMetrics (SM_CXSCREEN); //размер экрана
  		    	maxY = GetSystemMetrics (SM_CYSCREEN);
  			hdc = GetDC(hwnd); //получить контекст устройства
  			memdc1 = CreateCompatibleDC(hdc); //совместимый контекст в ОЗУ
  			hbit = CreateCompatibleBitmap (hdc,maxX,maxY); //место в памяти
									// для содержимого окна
  			SelectObject (memdc1, hbit); // выбрать bitmap
  			hbrush = (HBRUSH) GetStockObject (GRAY_BRUSH); //получить кисть
  			SelectObject (memdc1, hbrush); // выбрать кисть
  			PatBlt (memdc1, 0,0, maxX,maxY,PATCOPY); //закрасить кистью копию окна в ОЗУ 
  			ReleaseDC (hwnd, hdc); //освободить контекст        
				
			break;
  		
		case WM_PAINT: //отображение в окне на экране
			
			hdc = BeginPaint(hwnd,&paintstruct);
			BitBlt (hdc, 0, 0, maxX,maxY, memdc1,0,0,SRCCOPY);			
			EndPaint(hwnd,&paintstruct);
			break;

  		case WM_COMMAND:
  			switch (LOWORD(wParam))
  			{
  				
				case ID_PAINT:
					ShowWindow(hwndson,SW_SHOWNORMAL);
				break;	
				
				case ID_HIDE:
					ShowWindow(hwndson,SW_HIDE);
				break;	

				case ID_CLEAR:
					eraser=TRUE;
					InvalidateRect(hwndson,NULL,TRUE);
					UpdateWindow (hwndson);
				break;
				
				case ID_ABOUT:
  					MessageBox(hwnd, "Автор: Тигран Калайджян\n212 группа",
  					"Кратко о программе",MB_OK | MB_ICONINFORMATION);
  					break;
  			
				
  				case ID_EXIT:
  				    if (MessageBox(hwnd, "Вы уверены?",
  					"Выход",MB_YESNO | MB_ICONQUESTION) == IDYES)
  				    	DestroyWindow(hwnd);
  					break;
  			    
				case ID_INSDATA:

					hwndDlg=CreateDialog(GetModuleHandle(NULL),MAKEINTRESOURCE(IDD_TIGDG),hwnd,DialogProc);
					if (hwnd != NULL) ShowWindow(hwndDlg,SW_SHOW);
					else MessageBox(hwnd, "Ошибка создания диалога", "ERROR",MB_OK | MB_ICONINFORMATION);
  					
					break;
  			    
				case ID_CHECK:
				for (i=1;i<=dimmas-1;i++)
				for (j=1;j<=dimmas-1;j++)
				{
					matrix1[i][j]=matrix[i][j];
					x=i-(dimmas-1)/2;
					y=-j+(dimmas-1)/2;
					r=sqrt(x*x+y*y)/((dimmas-1)/2); // <- нормировка !!!
					if (x!=0) theta=atan2(y,x); else
					{
						if (y>0) theta=PI/2; else theta=3*PI/2; //предусмотреть вариант для ПИ
					}
					matrix[i][j]=0;
					for (s=0;s<number_contacts;s++) matrix[i][j]+=voltages(s+1)*response(s,r,theta);

					if (r>1) matrix[i][j]=0;
				}
				InvalidateRect(hwndson,NULL,TRUE);
				UpdateWindow (hwndson);
				MessageBox(hwnd, "Построение закончено", "OK",MB_OK | MB_ICONINFORMATION);
				break;
				
				case ID_KOMPENS:
				for (i=1;i<=dimmas-1;i++)
					for (j=1;j<=dimmas-1;j++)
						matrix[i][j]-=matrix1[i][j];
				eraser1=FALSE;
				/////////////////////////////////////////////
				//SetTextColor(hdc,255);
				//TextOut(hdc,100,100,mybuf,sizeof(mybuf));
				InvalidateRect(hwndson,NULL,TRUE);
				UpdateWindow (hwndson);
				InvalidateRect(hwnd,NULL,TRUE);
				UpdateWindow (hwnd);
				sprintf((LPTSTR)&mybuf,"среднеквадратичная ошибка компенсации\n %.5f длин волн\n\n \"сглаживание\" %.1f раз",
					findaver(),nachs/(findaver()));
				MessageBox(hwnd, mybuf, "Результат компенсации",MB_OK);
				
				break;

				case ID_APPROX:
					int i,j,k;
					column.setbounds(1, number_contacts);
					voltages.setbounds(1, number_contacts);
					square_matrix.setbounds(1,number_contacts,1,number_contacts);

					for (j=1;j<=number_contacts;j++)
					{
						column(j)=0;
						voltages(j)=0;
						for (theta=0;theta<2*PI;theta+=(2*PI/step))
							for (ro=0.0001;ro<1;ro+=1.0/step)// <- STEP!
							column(j)+=matrix[(int)((dimmas-1)*(1+ro*cos(theta))/2)][(int)
									((dimmas-1)*(1-ro*sin(theta))/2)]*response(j-1,ro,theta);
						
						for (i=1;i<=number_contacts;i++)
						{
							square_matrix(j,i)=0;
							for (theta=0;theta<2*PI;theta+=2*PI/step)
								for (ro=0.0001;ro<1;ro+=1.0/step) 
									square_matrix(j,i)+=response(i-1,ro,theta)*response(j-1,ro,theta);
						}
					}
					solvesystem(square_matrix,column,number_contacts,voltages);
					for (j=1;j<=number_contacts;j++)
						if (fabs(voltages(j))>300) voltages(j)=300*fabs(voltages(j))/voltages(j);
					sprintf((LPTSTR)&mybuf,"1-й контакт: %.2f\n2-й контакт: %.2f\n3-й контакт: %.2f\n4-й контакт: 
						%.2f\n5-й контакт: %.2f\n6-й контакт: %.2f\n7-й контакт: %.2f\n8-й контакт: %.2f\n9-й 
						контакт: %.2f\n10-й контакт: %.2f\n11-й контакт: %.2f\n12-й контакт: %.2f\n13-й контакт: 
						%.2f\n14-й контакт: %.2f\n15-й контакт: %.2f\n16-й контакт: %.2f\n17-й контакт: %.2f	\n”,voltages(1),voltages(2),voltages(3),voltages(4),voltages(5),voltages(6),voltages(7),voltages(8),voltages(9),voltages(10),voltages(11),
						voltages(12),voltages(13),voltages(14),voltages(15),voltages(16),voltages(17));
					MessageBox(hwnd, mybuf, "Напряжения (В)",MB_OK);
					hMenu=GetMenu (hwnd);
					EnableMenuItem(hMenu,ID_CHECK,MF_ENABLED);
					EnableMenuItem(hMenu,ID_KOMPENS,MF_ENABLED);
					break;
  			}
  			break;
  		            	
		case WM_TIMER:
			eraser=TRUE;
			InvalidateRect(hwndson,NULL,TRUE);
			UpdateWindow (hwndson);
				
			break;

  		case WM_DESTROY:
//  		    WritePrivateProfileString("main","global_a",&mybuf,iniad);
			PostQuitMessage(0);
  			break;
  		
  		default:
  			return DefWindowProc (hwnd, message, wParam, lParam);
  	}
  	return 0;
  }							
  	
BOOL CALLBACK DialogProc(HWND hDlg, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
	int i,j,v;
	int x,y;
	double r;
	long double theta,psi;
    	switch (uMsg)
  	{
  		case WM_INITDIALOG:

			 KillTimer(hwnd,NULL);
  			 sprintf((LPTSTR)&mybuf,"%.3f",zernike[1]);
			 SetDlgItemText(hDlg,IDC_EDIT1,(LPTSTR)&mybuf);
			 sprintf((LPTSTR)&mybuf,"%.3f",zernike[2]);
			 SetDlgItemText(hDlg,IDC_EDIT2,(LPTSTR)&mybuf);
			 sprintf((LPTSTR)&mybuf,"%.3f",zernike[3]);
			 SetDlgItemText(hDlg,IDC_EDIT3,(LPTSTR)&mybuf);
             		 sprintf((LPTSTR)&mybuf,"%.3f",zernike[4]);
			 SetDlgItemText(hDlg,IDC_EDIT4,(LPTSTR)&mybuf);
             		 sprintf((LPTSTR)&mybuf,"%.3f",zernike[5]);
			 SetDlgItemText(hDlg,IDC_EDIT5,(LPTSTR)&mybuf);
             		 sprintf((LPTSTR)&mybuf,"%.3f",zernike[6]);
			 SetDlgItemText(hDlg,IDC_EDIT6,(LPTSTR)&mybuf);
               		 sprintf((LPTSTR)&mybuf,"%.3f",zernike[7]);
			 SetDlgItemText(hDlg,IDC_EDIT7,(LPTSTR)&mybuf);
			 sprintf((LPTSTR)&mybuf,"%.3f",zernike[10]);
			 SetDlgItemText(hDlg,IDC_EDIT8,(LPTSTR)&mybuf);
			 SetFocus(GetDlgItem(hDlg,IDC_EDIT1));
			 break;

		case WM_CLOSE:
			 InvalidateRect(hwndson,NULL,TRUE);
  			 UpdateWindow (hwndson);
			 EndDialog(hDlg,NULL);
			 hwndDlg=0;
//			 SetTimer(hwnd,NULL,global_t,NULL);

			 break;

  		case WM_COMMAND:
  			if (HIWORD(wParam)==BN_CLICKED)
			{
				switch (wParam)
  				{
  				case IDC_CANCEL:
				SendMessage(hDlg,WM_CLOSE,NULL,NULL);
				break;

				case IDC_BUTTON1:
				srand((unsigned)time( NULL ));//выбираем семя рандомизации
				tmp1=(rand()%10);
				tmp2=(rand()%10);
				tmp3=1+(rand()%10);

				for (i=1;i<=dimmas-1;i++)
				for (j=1;j<=dimmas-1;j++)
				{
					x=i-(dimmas-1)/2;
					y=-j+(dimmas-1)/2;
					r=sqrt(x*x+y*y)/((dimmas-1)/2); // <- нормировка !!!
					if (x!=0) theta=atan2(y,x); else
					{
						if (y>0) theta=PI/2; else theta=3*PI/2; //предусмотреть вариант для ПИ
					}
					matrix[i][j]=0.01*tmp3*sin(tmp3*100/(r+0.0001));

					if (r>1) matrix[i][j]=0;
				}

				SendMessage(hDlg,WM_CLOSE,NULL,NULL);
				nachs=findaver();
				break;


				case IDOK:

				GetDlgItemText(hDlg,IDC_EDIT1,(LPTSTR)&mybuf,20);
				zernike[1] = atof((LPTSTR)&mybuf);
				GetDlgItemText(hDlg,IDC_EDIT2,(LPTSTR)&mybuf,20);
				zernike[2] = atof((LPTSTR)&mybuf);
				GetDlgItemText(hDlg,IDC_EDIT3,(LPTSTR)&mybuf,20);
				zernike[3] = atof((LPTSTR)&mybuf);
				GetDlgItemText(hDlg,IDC_EDIT4,(LPTSTR)&mybuf,20);
				zernike[4] = atof((LPTSTR)&mybuf);
				GetDlgItemText(hDlg,IDC_EDIT5,(LPTSTR)&mybuf,20);
				zernike[5] = atof((LPTSTR)&mybuf);
				GetDlgItemText(hDlg,IDC_EDIT6,(LPTSTR)&mybuf,20);
				zernike[6] = atof((LPTSTR)&mybuf);
				GetDlgItemText(hDlg,IDC_EDIT7,(LPTSTR)&mybuf,20);
				zernike[7] = atof((LPTSTR)&mybuf);
				GetDlgItemText(hDlg,IDC_EDIT8,(LPTSTR)&mybuf,20);
				zernike[10] = atof((LPTSTR)&mybuf);

				// записываем значения в массив
				//for(i=0;i<=6;i++) norma+= (zernike[i]);
				//norma=(zernike[0]*2+zernike[1]*2+zernike[2]*sqrt(6)+zernike[3]*sqrt(6)+
						zernike[4]*2*sqrt(3)+sqrt(8)*zernike[5]+sqrt(8)*zernike[6]);
				for (i=1;i<=dimmas-1;i++)
					for (j=1;j<=dimmas-1;j++)
					{
						x=i-(dimmas-1)/2;
						y=-j+(dimmas-1)/2;
						r=sqrt(x*x+y*y)/((dimmas-1)/2); // <- нормировка !!!
						if (x!=0) theta=atan2(y,x); else
						{
							if (y>0) theta=PI/2; else theta=3*PI/2; //предусмотреть вариант для ПИ
						}
						matrix[i][j]=0;
						for (v=0;v<16;v++) matrix[i][j]+=zernike[v]*polynomial(v+1,r,theta);
						if (r>1) matrix[i][j]=0;
					}

				///////////////////////////////
				SendMessage(hDlg,WM_CLOSE,NULL,NULL);
/*				if (global_a != 0) 
				{
					if ((global_a > 1.5)||(global_a < 0.5)) 
					{
						global_a=old;
						MessageBox(hwnd, "Выбранное значение не лежит в промежутке [0.5; 1.5]
							\nвведите что-нибудь другое.", "ERROR",MB_OK | MB_ICONERROR);
					}
					else SendMessage(hDlg,WM_CLOSE,NULL,NULL);
				}
				else 
				{
					global_a=old;
					MessageBox(hwnd, "Выбранное значение - ноль или не число,\nвведите что-нибудь другое.", "ERROR",MB_OK | MB_ICONERROR);
				}
*/
				nachs=findaver();
				break;
			}
		}
  			        		
  		default:
  		return 0;
	}
return 1;  
}							

LRESULT CALLBACK WindowSon(HWND hwndson,UINT message,
  							WPARAM wParam,LPARAM lParam)
  {
  	HDC hdc;
  	HBITMAP hbit;
  	HBRUSH hbrush;
  	int x,y;

  	static int maxX, maxY;  	
  	
  	switch (message)
  	{
  		case WM_CREATE: // событие, возникающее при создании окна

		    maxX = GetSystemMetrics (SM_CXSCREEN); //размер экрана
  		    maxY = GetSystemMetrics (SM_CYSCREEN);
  			hdc = GetDC(hwndson); //получить контекст устройства
  			memdc = CreateCompatibleDC(hdc); //совместимый контекст в ОЗУ
  			hbit = CreateCompatibleBitmap (hdc,maxX,maxY); //место в памяти
									// для содержимого окна
  			SelectObject (memdc, hbit); // выбрать bitmap
  			hbrush = (HBRUSH) GetStockObject (GRAY_BRUSH); //получить кисть
  			SelectObject (memdc, hbrush); // выбрать кисть
  			PatBlt (memdc, 0,0, maxX,maxY,PATCOPY); //закрасить кистью копию окна в ОЗУ 
  			ReleaseDC (hwndson, hdc); //освободить контекст        
		
			
			break;
  		
  		case WM_PAINT: //отображение в окне на экране
  			
			if (!eraser)
			{
				norma=findnorm();				
				int prob;
				for (x=1;x<=dimmas-1;x++)
				for (y=1;y<=dimmas-1;y++)				
				{
					prob=(int)(128*(1+matrix[x][y]/norma));  // <- нормировка
					SetPixel(memdc, x,y,RGB(prob,prob,prob)); 
				}
			}
			else 
			{
				SelectObject (memdc, hbrush); // выбрать кисть
  			    PatBlt (memdc, 0,0, maxX,maxY,PATCOPY); //закрасить кистью копию окна в ОЗУ 
				eraser=FALSE;
  			}

			hdc = BeginPaint(hwndson,&paintstruct);
			BitBlt (hdc, 0, 0, maxX,maxY, memdc,0,0,SRCCOPY);			
			EndPaint(hwndson,&paintstruct);
			break;
	
  		default:
  			return DefWindowProc (hwndson, message, wParam, lParam);
  	}
  	return 0;
  }							

double polynomial(int n, double r, double theta)
{
	switch (n)
	{
		case 1: return 1; break;
		case 2: return 2*r*cos(theta); break;
		case 3: return 2*r*sin(theta); break;
		case 4: return sqrt(3)*(2*r*r-1); break;
		case 5: return sqrt(6)*r*r*sin(2*theta); break;
		case 6: return sqrt(6)*r*r*cos(2*theta); break;
		case 7: return sqrt(8)*(3*r*r*r-2*r)*sin(theta); break;
		case 8: return sqrt(8)*(3*r*r*r-2*r)*cos(theta); break;
		case 9: return sqrt(8)*pow(r,3)*sin(3*theta); break;
		case 10: return sqrt(8)*pow(r,3)*cos(3*theta); break;
		case 11: return sqrt(5)*(6*pow(r,4)-6*r*r+1); break;
		case 12: return sqrt(10)*(4*pow(r,4)-3*r*r)*cos(2*theta); break;
		case 13: return sqrt(10)*(4*pow(r,4)-3*r*r)*sin(2*theta); break;
		case 14: return sqrt(10)*pow(r,4)*cos(4*theta); break;
		case 15: return sqrt(10)*pow(r,4)*sin(4*theta); break;
		case 16: return sqrt(12)*(10*pow(r,5)-12*pow(r,3)+3*r)*cos(theta); break;
	}
}

/* функции отклика. контакты нумеруются с НУЛЯ ! */
long double response(int kontakt, double r, double theta)
{
	long double result=0;
	int kou;
	for(kou=0;kou<=12;kou++) result+=response[kontakt][kou]*polynomial(kou+4,r,theta);
	return result;
}

/*
Решение системы  линейных  уравнений  с  матрицей  системы,  заданной  LU-
разложением.

Алгоритм решает систему линейных уравнений,  матрица  которой  задана  LU-
разложением. В случае, если  в  процессе  решения  произойдет  деление  на
ноль, возвращается значение  False,  обозначающее,  что  система  является
вырожденной.

Алгоритм решает только системы уравнений с квадратной матрицей.

Входные параметры:
    A   -   LU-разложение матрицы системы в упакованной  форме  (результат
	    работы подпрограммы LUDecomposition).
    Pivots- таблица  перестановок  строк  (результат  работы  подпрограммы
	    LUDecomposition).
    B   -   правая часть системы. Массив с нумерацией элементов [1..N]
    N   -   размерность системы.

Выходные параметры:
    X   -   решение системы. Массив с нумерацией элементов [1..N]

Результат:
    True, если система не вырождена (но, возможно, близка к вырожденной).
    False, если система вырождена. В таком случае X не содержит решение.
*************************************************************************/
bool solvesystemlu(const ap::real_2d_array& a,
     const ap::integer_1d_array& pivots,
     ap::real_1d_array b,
     int n,
     ap::real_1d_array& x)
{
    bool result;
    ap::real_1d_array y;
    int i;
    int j;
    double v;
    int ip1;
    int im1;

    y.setbounds(1, n);
    x.setbounds(1, n);
    result = true;
    for(i = 1; i <= n; i++)
    {
        if( a(i,i)==0 )
        {
            result = false;
            return result;
        }
    }
    for(i = 1; i <= n; i++)
    {
        if( pivots(i)!=i )
        {
            v = b(i);
            b(i) = b(pivots(i));
            b(pivots(i)) = v;
        }
    }
    y(1) = b(1);
    for(i = 2; i <= n; i++)
    {
        im1 = i-1;
        v = ap::vdotproduct(a.getrow(i, 1, im1), y.getvector(1, im1));
        y(i) = b(i)-v;
    }
    x(n) = y(n)/a(n,n);
    for(i = n-1; i >= 1; i--)
    {
        ip1 = i+1;
        v = ap::vdotproduct(a.getrow(i, ip1, n), x.getvector(ip1, n));
        x(i) = (y(i)-v)/a(i,i);
    }
    return result;
}

/*************************************************************************
Решение системы  линейных  уравнений

Алгоритм решает систему линейных уравнений с использованием LU-разложения.
Алгоритм решает только системы уравнений с квадратной матрицей.

Входные параметры:
    A   -   Матрица системы.
            Массив с нумерацией элементов [1..N, 1..N].
    B   -   Правая часть.
            Массив с нумерацией элементов [1..N, 1..N].
    N   -   размерность системы.

Выходные параметры:
    X   -   решение системы. Массив с нумерацией элементов [1..N]

Результат:
    True, если система не вырождена (но, возможно, близка к вырожденной).
    False, если система вырождена. В таком случае X не содержит решение.
*************************************************************************/
bool solvesystem(ap::real_2d_array a,
     ap::real_1d_array b,
     int n,
     ap::real_1d_array& x)
{
    bool result;
    ap::integer_1d_array pivots;
    int i;

    ludecomposition(a, n, n, pivots);
    result = solvesystemlu(a, pivots, b, n, x);
    return result;
}

/*************************************************************************
LU-разложение матрицы общего вида размера M x N

Подпрограмма вычисляет LU-разложение прямоугольной матрицы общего  вида  с
частичным выбором ведущего элемента (с перестановками строк).

Входные параметры:
    A       -   матрица A. Нумерация элементов: [1..M, 1..N]
    M       -   число строк в матрице A
    N       -   число столбцов в матрице A

Выходные параметры:
    A       -   матрицы L и U в компактной форме (см. ниже).
		Нумерация элементов: [1..M, 1..N]
    Pivots  -   матрица перестановок в компактной форме (см. ниже).
		Нумерация элементов: [1..Min(M,N)]

Матрица A представляется, как A = P * L * U, где P - матрица перестановок,
матрица L - нижнетреугольная (или нижнетрапецоидальная, если M>N) матрица,
U - верхнетреугольная (или верхнетрапецоидальная, если M<N) матрица.

Рассмотрим разложение более подробно на примере при M=4, N=3:

				   (  1          )    ( U11 U12 U13  )
A = P1 * P2 * P3 * ( L21  1      )  * (     U22 U23  )
				   ( L31 L32  1  )    (         U33  )
		           ( L41 L42 L43 )

Здесь матрица L  имеет  размер  M  x  Min(M,N),  матрица  U  имеет  размер
Min(M,N) x N, матрица  P(i)  получается  путем  перестановки  в  единичной
матрице размером M x M строк с номерами I и Pivots[I]

Результатом работы алгоритма являются массив Pivots  и  следующая матрица,
замещающая  матрицу  A,  и  сохраняющая  в компактной форме матрицы L и U
(пример приведен для M=4, N=3):

 ( U11 U12 U13 )
 ( L21 U22 U23 )
 ( L31 L32 U33 )
 ( L41 L42 L43 )

Как видно, единичная диагональ матрицы L  не  сохраняется.
Если N>M, то соответственно меняются размеры матриц и расположение
элементов.
*************************************************************************/
void ludecomposition(ap::real_2d_array& a,
     int m,
     int n,
     ap::integer_1d_array& pivots)
{
    int i;
    int j;
    int jp;
    ap::real_1d_array t1;
    double s;

    pivots.setbounds(1, ap::minint(m, n));
    t1.setbounds(1, ap::maxint(m, n));
    ap::ap_error::make_assertion(m>=0&&n>=0);
    if( m==0||n==0 )
    {
        return;
    }
    for(j = 1; j <= ap::minint(m, n); j++)
    {
        jp = j;
        for(i = j+1; i <= m; i++)
        {
            if( fabs(a(i,j))>fabs(a(jp,j)) )
            {
                jp = i;
            }
        }
        pivots(j) = jp;
        if( a(jp,j)!=0 )
        {
            if( jp!=j )
            {
                ap::vmove(t1.getvector(1, n), a.getrow(j, 1, n));
                ap::vmove(a.getrow(j, 1, n), a.getrow(jp, 1, n));
                ap::vmove(a.getrow(jp, 1, n), t1.getvector(1, n));
            }
            if( j<m )
            {
                jp = j+1;
                s = double(1)/double(a(j,j));
                ap::vmul(a.getcolumn(j, jp, m), s);
            }
        }
        if( j<ap::minint(m, n) )
        {
            jp = j+1;
            for(i = j+1; i <= m; i++)
            {
                s = a(i,j);
                ap::vsub(a.getrow(i, jp, n), a.getrow(j, jp, n), s);
            }
        }
    }
}

long double findnorm(void)
{
long double max=0;
int i,j;
for (i=0;i<dimmas;i++)
	for (j=0;j<dimmas;j++)
		if (fabs(matrix[i][j])>max) max=fabs(matrix[i][j]);
		return max;
}

long double findaver(void)
{
const step1=100;
long double average=0, sdev=0,theta,ro;
for (theta=0;theta<2*PI;theta+=(2*PI/step1))
		for (ro=0.0001;ro<1;ro+=1.0/step1)
		average+=matrix[(int)((dimmas-1)*(1+ro*cos(theta))/2)][(int)((dimmas-1)*(1-ro*sin(theta))/2)]/pow(step1,2);
for (theta=0;theta<2*PI;theta+=(2*PI/step1))
		for (ro=0.0001;ro<1;ro+=1.0/step1)
sdev+=pow(average - matrix[(int)((dimmas-1)*(1+ro*cos(theta))/2)][(int)((dimmas-1)*(1-ro*sin(theta))/2)],2)/pow(step1,4);
sdev=pow(sdev,0.5);
		return sdev;
}
