#ifndef AP_H
#define AP_H

#include <stdlib.h>
#include <math.h>

namespace ap
{

/********************************************************************
КЛАСС ИСКЛЮЧЕНИЯ, КОТОРОЕ ВЫБРАСЫВАЕТСЯ ПРИ РАЗЛИЧНЫХ ОШИБКАХ.

Текущая версия не позволяет определить причину, по которой исключение
было сгенерировано.
********************************************************************/
class ap_error
{
public:
    static void make_assertion(bool bClause)
        { if(!bClause) throw ap_error(); };
private:
};


/********************************************************************
ОПЕРАЦИЯ ЛОГИЧЕСКОЕ ИСКЛЮЧАЮЩЕЕ ИЛИ
********************************************************************/
static inline bool XOR(bool lhs, bool rhs)
{
    return (lhs && !rhs) || ( !lhs && rhs);
}


/********************************************************************
Класс-шаблон, определяющий вектор в памяти.
Используется базовыми подпрограммами линейной алгебры.

Вектор состоит из Length элементов типа T, расположенных с шагом Step
начиная с элемента, заданного указателем Data.

Этот класс поддерживает только константный доступ к вектору, запрещая
изменять его. 
********************************************************************/
template<class T>
class const_raw_vector
{
public:
    /****************************************************************
    создание вектора
    ****************************************************************/
    const_raw_vector(const T *Data, int Length, int Step):
        pData(const_cast<T*>(Data)),iLength(Length),iStep(Step){};

    /****************************************************************
    получение константного указателя на вектор
    ****************************************************************/
    const T* GetData() const
    { return pData; };

    /****************************************************************
    получение числа компонент вектора
    ****************************************************************/
    int GetLength() const
    { return iLength; };

    /****************************************************************
    получение шага между компонентами вектора
    ****************************************************************/
    int GetStep() const
    { return iStep; };
protected:
    T       *pData;
    int     iLength, iStep;
};


/********************************************************************
Класс-шаблон, определяющий вектор в памяти.
Наследует от const_raw_vector.
Используется базовыми подпрограммами линейной алгебры.

Вектор состоит из Length элементов типа T, расположенных с шагом Step
начиная с элемента, заданного указателем Data.

Этот класс поддерживает произвольный доступ к вектору,  разрешая  как
чтение, так и изменение компонент вектора.
********************************************************************/
template<class T>
class raw_vector : public const_raw_vector<T>
{
public:
    /****************************************************************
    создание вектора
    ****************************************************************/
    raw_vector(T *Data, int Length, int Step):const_raw_vector<T>(Data, Length, Step){};

    /****************************************************************
    получение указателя на вектор
    ****************************************************************/
    T* GetData()
    { return pData; };
};


/********************************************************************
Скалярное произведение векторов

Функция возвращает скалярное произведение векторов v1 и v2. Для
вычисления используются развернутые циклы.
********************************************************************/
template<class T>
T vdotproduct(const_raw_vector<T> v1, const_raw_vector<T> v2)
{
    ap_error::make_assertion(v1.GetLength()==v2.GetLength());
    if( v1.GetStep()==1 && v2.GetStep()==1 )
    {
        //
        // fast
        //
        T r = 0;
        const T *p1 = v1.GetData();
        const T *p2 = v2.GetData();
        int imax = v1.GetLength()/4;
        int i;
        for(i=imax; i!=0; i--)
        {
            r += (*p1)*(*p2) + p1[1]*p2[1] + p1[2]*p2[2] + p1[3]*p2[3];
            p1+=4;
            p2+=4;
        }
        for(i=0; i<v1.GetLength()%4; i++)
            r += (*(p1++))*(*(p2++));
        return r;
    }
    else
    {
        //
        // general
        //
        int offset11 = v1.GetStep(), offset12 = 2*offset11, offset13 = 3*offset11, offset14 = 4*offset11;
        int offset21 = v2.GetStep(), offset22 = 2*offset21, offset23 = 3*offset21, offset24 = 4*offset21;
        T r = 0;
        const T *p1 = v1.GetData();
        const T *p2 = v2.GetData();
        int imax = v1.GetLength()/4;
        int i;
        for(i=0; i<imax; i++)
        {
            r += (*p1)*(*p2) + p1[offset11]*p2[offset21] + p1[offset12]*p2[offset22] + p1[offset13]*p2[offset23];
            p1+=offset14;
            p2+=offset24;
        }
        for(i=0; i<v1.GetLength()%4; i++)
        {
            r += (*p1)*(*p2);
            p1+=offset11;
            p2+=offset21;
        }
        return r;
    }
}


/********************************************************************
Копирование вектора

Функция копирует вектор vsrc в вектор vdst.
********************************************************************/
template<class T>
void vmove(raw_vector<T> vdst, const_raw_vector<T> vsrc)
{
    ap_error::make_assertion(vdst.GetLength()==vsrc.GetLength());
    if( vdst.GetStep()==1 && vsrc.GetStep()==1 )
    {
        //
        // fast
        //
        T *p1 = vdst.GetData();
        const T *p2 = vsrc.GetData();
        int imax = vdst.GetLength()/2;
        int i;
        for(i=imax; i!=0; i--)
        {
            *p1 = *p2;
            p1[1] = p2[1];
            p1 += 2;
            p2 += 2;
        }
        if(vdst.GetLength()%2 != 0)
            *p1 = *p2;
        return;
    }
    else
    {
        //
        // general
        //
        int offset11 = vdst.GetStep(), offset12 = 2*offset11, offset13 = 3*offset11, offset14 = 4*offset11;
        int offset21 = vsrc.GetStep(), offset22 = 2*offset21, offset23 = 3*offset21, offset24 = 4*offset21;
        T *p1 = vdst.GetData();
        const T *p2 = vsrc.GetData();
        int imax = vdst.GetLength()/4;
        int i;
        for(i=0; i<imax; i++)
        {
            *p1 = *p2;
            p1[offset11] = p2[offset21];
            p1[offset12] = p2[offset22];
            p1[offset13] = p2[offset23];
            p1 += offset14;
            p2 += offset24;
        }
        for(i=0; i<vdst.GetLength()%4; i++)
        {
            *p1 = *p2;
            p1 += vdst.GetStep();
            p2 += vsrc.GetStep();
        }
        return;
    }
}


/********************************************************************
Копирование вектора с умножением на -1

Функция копирует вектор vsrc в вектор vdst, умножая каждую его
компоненту на -1.
********************************************************************/
template<class T>
void vmoveneg(raw_vector<T> vdst, const_raw_vector<T> vsrc)
{
    ap_error::make_assertion(vdst.GetLength()==vsrc.GetLength());
    if( vdst.GetStep()==1 && vsrc.GetStep()==1 )
    {
        //
        // fast
        //
        T *p1 = vdst.GetData();
        const T *p2 = vsrc.GetData();
        int imax = vdst.GetLength()/2;
        int i;
        for(i=0; i<imax; i++)
        {
            *p1 = -*p2;
            p1[1] = -p2[1];
            p1 += 2;
            p2 += 2;
        }
        if(vdst.GetLength()%2 != 0)
            *p1 = -*p2;
        return;
    }
    else
    {
        //
        // general
        //
        int offset11 = vdst.GetStep(), offset12 = 2*offset11, offset13 = 3*offset11, offset14 = 4*offset11;
        int offset21 = vsrc.GetStep(), offset22 = 2*offset21, offset23 = 3*offset21, offset24 = 4*offset21;
        T *p1 = vdst.GetData();
        const T *p2 = vsrc.GetData();
        int imax = vdst.GetLength()/4;
        int i;
        for(i=imax; i!=0; i--)
        {
            *p1 = -*p2;
            p1[offset11] = -p2[offset21];
            p1[offset12] = -p2[offset22];
            p1[offset13] = -p2[offset23];
            p1 += offset14;
            p2 += offset24;
        }
        for(i=0; i<vdst.GetLength()%4; i++)
        {
            *p1 = -*p2;
            p1 += vdst.GetStep();
            p2 += vsrc.GetStep();
        }
        return;
    }
}


/********************************************************************
Копирование вектора с умножением на произвольное число

Функция копирует вектор vsrc в вектор vdst, умножая каждую компоненту
вектора vsrc на alpha.
********************************************************************/
template<class T, class T2>
void vmove(raw_vector<T> vdst, const_raw_vector<T> vsrc, T2 alpha)
{
    ap_error::make_assertion(vdst.GetLength()==vsrc.GetLength());
    if( vdst.GetStep()==1 && vsrc.GetStep()==1 )
    {
        //
        // fast
        //
        T *p1 = vdst.GetData();
        const T *p2 = vsrc.GetData();
        int imax = vdst.GetLength()/4;
        int i;
        for(i=imax; i!=0; i--)
        {
            *p1 = alpha*(*p2);
            p1[1] = alpha*p2[1];
            p1[2] = alpha*p2[2];
            p1[3] = alpha*p2[3];
            p1 += 4;
            p2 += 4;
        }
        for(i=0; i<vdst.GetLength()%4; i++)
            *(p1++) = alpha*(*(p2++));
        return;
    }
    else
    {
        //
        // general
        //
        int offset11 = vdst.GetStep(), offset12 = 2*offset11, offset13 = 3*offset11, offset14 = 4*offset11;
        int offset21 = vsrc.GetStep(), offset22 = 2*offset21, offset23 = 3*offset21, offset24 = 4*offset21;
        T *p1 = vdst.GetData();
        const T *p2 = vsrc.GetData();
        int imax = vdst.GetLength()/4;
        int i;
        for(i=0; i<imax; i++)
        {
            *p1 = alpha*(*p2);
            p1[offset11] = alpha*p2[offset21];
            p1[offset12] = alpha*p2[offset22];
            p1[offset13] = alpha*p2[offset23];
            p1 += offset14;
            p2 += offset24;
        }
        for(i=0; i<vdst.GetLength()%4; i++)
        {
            *p1 = alpha*(*p2);
            p1 += vdst.GetStep();
            p2 += vsrc.GetStep();
        }
        return;
    }
}


/********************************************************************
Добавление вектора

Функция добавляет вектор vsrc к вектору vdst
********************************************************************/
template<class T>
void vadd(raw_vector<T> vdst, const_raw_vector<T> vsrc)
{
    ap_error::make_assertion(vdst.GetLength()==vsrc.GetLength());
    if( vdst.GetStep()==1 && vsrc.GetStep()==1 )
    {
        //
        // fast
        //
        T *p1 = vdst.GetData();
        const T *p2 = vsrc.GetData();
        int imax = vdst.GetLength()/4;
        int i;
        for(i=imax; i!=0; i--)
        {
            *p1 += *p2;
            p1[1] += p2[1];
            p1[2] += p2[2];
            p1[3] += p2[3];
            p1 += 4;
            p2 += 4;
        }
        for(i=0; i<vdst.GetLength()%4; i++)
            *(p1++) += *(p2++);
        return;
    }
    else
    {
        //
        // general
        //
        int offset11 = vdst.GetStep(), offset12 = 2*offset11, offset13 = 3*offset11, offset14 = 4*offset11;
        int offset21 = vsrc.GetStep(), offset22 = 2*offset21, offset23 = 3*offset21, offset24 = 4*offset21;
        T *p1 = vdst.GetData();
        const T *p2 = vsrc.GetData();
        int imax = vdst.GetLength()/4;
        int i;
        for(i=0; i<imax; i++)
        {
            *p1 += *p2;
            p1[offset11] += p2[offset21];
            p1[offset12] += p2[offset22];
            p1[offset13] += p2[offset23];
            p1 += offset14;
            p2 += offset24;
        }
        for(i=0; i<vdst.GetLength()%4; i++)
        {
            *p1 += *p2;
            p1 += vdst.GetStep();
            p2 += vsrc.GetStep();
        }
        return;
    }
}


/********************************************************************
Добавление вектора с умножением на произвольное число

Функция добавляет вектор vsrc к вектору vdst, умножая каждую компоненту
вектора vsrc на alpha
********************************************************************/
template<class T, class T2>
void vadd(raw_vector<T> vdst, const_raw_vector<T> vsrc, T2 alpha)
{
    ap_error::make_assertion(vdst.GetLength()==vsrc.GetLength());
    if( vdst.GetStep()==1 && vsrc.GetStep()==1 )
    {
        //
        // fast
        //
        T *p1 = vdst.GetData();
        const T *p2 = vsrc.GetData();
        int imax = vdst.GetLength()/4;
        int i;
        for(i=imax; i!=0; i--)
        {
            *p1 += alpha*(*p2);
            p1[1] += alpha*p2[1];
            p1[2] += alpha*p2[2];
            p1[3] += alpha*p2[3];
            p1 += 4;
            p2 += 4;
        }
        for(i=0; i<vdst.GetLength()%4; i++)
            *(p1++) += alpha*(*(p2++));
        return;
    }
    else
    {
        //
        // general
        //
        int offset11 = vdst.GetStep(), offset12 = 2*offset11, offset13 = 3*offset11, offset14 = 4*offset11;
        int offset21 = vsrc.GetStep(), offset22 = 2*offset21, offset23 = 3*offset21, offset24 = 4*offset21;
        T *p1 = vdst.GetData();
        const T *p2 = vsrc.GetData();
        int imax = vdst.GetLength()/4;
        int i;
        for(i=0; i<imax; i++)
        {
            *p1 += alpha*(*p2);
            p1[offset11] += alpha*p2[offset21];
            p1[offset12] += alpha*p2[offset22];
            p1[offset13] += alpha*p2[offset23];
            p1 += offset14;
            p2 += offset24;
        }
        for(i=0; i<vdst.GetLength()%4; i++)
        {
            *p1 += alpha*(*p2);
            p1 += vdst.GetStep();
            p2 += vsrc.GetStep();
        }
        return;
    }
}


/********************************************************************
Вычитание вектора

Функция вычитает вектор vsrc из вектора vdst
********************************************************************/
template<class T>
void vsub(raw_vector<T> vdst, const_raw_vector<T> vsrc)
{
    ap_error::make_assertion(vdst.GetLength()==vsrc.GetLength());
    if( vdst.GetStep()==1 && vsrc.GetStep()==1 )
    {
        //
        // fast
        //
        T *p1 = vdst.GetData();
        const T *p2 = vsrc.GetData();
        int imax = vdst.GetLength()/4;
        int i;
        for(i=imax; i!=0; i--)
        {
            *p1 -= *p2;
            p1[1] -= p2[1];
            p1[2] -= p2[2];
            p1[3] -= p2[3];
            p1 += 4;
            p2 += 4;
        }
        for(i=0; i<vdst.GetLength()%4; i++)
            *(p1++) -= *(p2++);
        return;
    }
    else
    {
        //
        // general
        //
        int offset11 = vdst.GetStep(), offset12 = 2*offset11, offset13 = 3*offset11, offset14 = 4*offset11;
        int offset21 = vsrc.GetStep(), offset22 = 2*offset21, offset23 = 3*offset21, offset24 = 4*offset21;
        T *p1 = vdst.GetData();
        const T *p2 = vsrc.GetData();
        int imax = vdst.GetLength()/4;
        int i;
        for(i=0; i<imax; i++)
        {
            *p1 -= *p2;
            p1[offset11] -= p2[offset21];
            p1[offset12] -= p2[offset22];
            p1[offset13] -= p2[offset23];
            p1 += offset14;
            p2 += offset24;
        }
        for(i=0; i<vdst.GetLength()%4; i++)
        {
            *p1 -= *p2;
            p1 += vdst.GetStep();
            p2 += vsrc.GetStep();
        }
        return;
    }
}


/********************************************************************
Вычитание вектора с умножением на произвольное число

Функция вычитает вектор vsrc из вектора vdst, умножая каждую компоненту
вектора vsrc на alpha
********************************************************************/
template<class T, class T2>
void vsub(raw_vector<T> vdst, const_raw_vector<T> vsrc, T2 alpha)
{
    vadd(vdst, vsrc, -alpha);
}


/********************************************************************
Умножение вектора

Функция умножает каждую компоненту вектора vdst на число alpha
********************************************************************/
template<class T, class T2>
void vmul(raw_vector<T> vdst, T2 alpha)
{
    if( vdst.GetStep()==1 )
    {
        //
        // fast
        //
        T *p1 = vdst.GetData();
        int imax = vdst.GetLength()/4;
        int i;
        for(i=imax; i!=0; i--)
        {
            *p1 *= alpha;
            p1[1] *= alpha;
            p1[2] *= alpha;
            p1[3] *= alpha;
            p1 += 4;
        }
        for(i=0; i<vdst.GetLength()%4; i++)
            *(p1++) *= alpha;
        return;
    }
    else
    {
        //
        // general
        //
        int offset11 = vdst.GetStep(), offset12 = 2*offset11, offset13 = 3*offset11, offset14 = 4*offset11;
        T *p1 = vdst.GetData();
        int imax = vdst.GetLength()/4;
        int i;
        for(i=0; i<imax; i++)
        {
            *p1 *= alpha;
            p1[offset11] *= alpha;
            p1[offset12] *= alpha;
            p1[offset13] *= alpha;
            p1 += offset14;
        }
        for(i=0; i<vdst.GetLength()%4; i++)
        {
            *p1 *= alpha;
            p1 += vdst.GetStep();
        }
        return;
    }
}


/********************************************************************
КЛАСС-ШАБЛОН ОДНОМЕРНОГО ДИНАМИЧЕСКОГО МАССИВА


ОПИСАНИЕ ЧЛЕНОВ КЛАССА:

---------------------------------------------------------------------
template_1d_array()

Создание пустого массива.

---------------------------------------------------------------------
~template_1d_array()

Удаление массива. При этом освобождается выделенная под массив память

---------------------------------------------------------------------
template_1d_array(const template_1d_array &rhs)

Создание копии массива. При этом выделяется отдельная область памяти,
в которую копируется содержимое массива-источника

---------------------------------------------------------------------
const template_1d_array& operator=(const template_1d_array &rhs)

Присваивание массива. При этом содержимое массива-приемника удаляется
и освобождается  выделенная под него память,  затем заново выделяется 
отдельная область памяти, в которую копируется содержимое источника.

---------------------------------------------------------------------
      T& operator()(int i)
const T& operator()(int i) const

Обращение к элементу массива с номером i

---------------------------------------------------------------------
void setbounds( int iLow, int iHigh )

Выделение памяти  под  массив.  При  этом  старое  содержимое массива
удаляется  и освобождается  выделенная под него память,  затем заново 
выделяется отдельная область памяти размера iHigh-iLow+1 элементов.

Нумерация элементов в новом массива начинается с iLow и заканчивается
iHigh. Содержимое нового массива не определено.

---------------------------------------------------------------------
void setcontent( int iLow, int iHigh, const T *pContent )

Метод  аналогичен  методу  setbounds()  за тем исключением, что после 
выделения памяти в неё копируется содержимое массива pContent[].

---------------------------------------------------------------------
      T* getcontent()
const T* getcontent() const

Метод позволяет получить указатель на содержимое массива. Данные,  на
которые указывает возвращенный указатель, можно изменять, и при  этом
изменится содержимое массива.

---------------------------------------------------------------------
int getlowbound()
int gethighbound()

Методы используются для  получения  информации  о  нижней  и  верхней
границах массива.


---------------------------------------------------------------------
raw_vector<T> getvector(int iStart, int iEnd)
const_raw_vector<T> getvector(int iStart, int iEnd) const

Методы используются базовыми подпрограммами линейной алгебры для
получения доступа к внутренней памяти массива. Методы возвращают
объект, содержащий в себе указатель на часть вектора (начиная
с элемента с индексом iStart и заканчивая индексом iEnd).

Если iEnd<iStart, то считается, что задан пустой вектор.
********************************************************************/
template<class T>
class template_1d_array
{
public:
    template_1d_array()
    {
        m_Vec=0;
        m_iVecSize = 0;
    };

    ~template_1d_array()
    {
        if(m_Vec)
            delete[] m_Vec;
    };

    template_1d_array(const template_1d_array &rhs)
    {
        m_iVecSize = rhs.m_iVecSize;
        m_iLow = rhs.m_iLow;
        m_iHigh = rhs.m_iHigh;
        if(rhs.m_Vec)
        {
            m_Vec = new T[m_iVecSize];
            #ifndef UNSAFE_MEM_COPY
            for(int i=0; i<m_iVecSize; i++)
                m_Vec[i] = rhs.m_Vec[i];
            #else
            memcpy(m_Vec, rhs.m_Vec, m_iVecSize*sizeof(T));
            #endif
        }
        else
            m_Vec=0;
    };

    const template_1d_array& operator=(const template_1d_array &rhs)
    {
        if( this==&rhs )
            return *this;

        m_iLow = rhs.m_iLow;
        m_iHigh = rhs.m_iHigh;
        m_iVecSize = rhs.m_iVecSize;
        if(m_Vec)
            delete[] m_Vec;
        if(rhs.m_Vec)
        {
            m_Vec = new T[m_iVecSize];
            #ifndef UNSAFE_MEM_COPY
            for(int i=0; i<m_iVecSize; i++)
                m_Vec[i] = rhs.m_Vec[i];
            #else
            memcpy(m_Vec, rhs.m_Vec, m_iVecSize*sizeof(T));
            #endif
        }
        else
            m_Vec=0;
        return *this;
    };

    const T& operator()(int i) const
    {
        #ifndef NO_AP_ASSERT
        ap_error::make_assertion(i>=m_iLow && i<=m_iHigh);
        #endif
        return m_Vec[ i-m_iLow ];
    };

    T& operator()(int i)
    {
        #ifndef NO_AP_ASSERT
        ap_error::make_assertion(i>=m_iLow && i<=m_iHigh);
        #endif
        return m_Vec[ i-m_iLow ];
    };

    void setbounds( int iLow, int iHigh )
    {
        if(m_Vec)
            delete[] m_Vec;
        m_iLow = iLow;
        m_iHigh = iHigh;
        m_iVecSize = iHigh-iLow+1;
        m_Vec = new T[m_iVecSize];
    };

    void setcontent( int iLow, int iHigh, const T *pContent )
    {
        setbounds(iLow, iHigh);
        for(int i=iLow; i<=iHigh; i++)
            (*this)(i) = pContent[i-iLow];
    };

    T* getcontent()
    {
        return m_Vec;
    };

    const T* getcontent() const
    {
        return m_Vec;
    };

    int getlowbound(int iBoundNum = 0) const
    {
        return m_iLow;
    };

    int gethighbound(int iBoundNum = 0) const
    {
        return m_iHigh;
    };

    raw_vector<T> getvector(int iStart, int iEnd)
    {
        #ifndef NO_AP_ASSERT
        if( iEnd>=iStart )
            ap_error::make_assertion(m_iLow<=iStart && iEnd<=m_iHigh);
        #endif
        if( iStart>iEnd )
            return raw_vector<T>(0, 0, 1);
        else
            return raw_vector<T>(m_Vec+iStart-m_iLow, iEnd-iStart+1, 1);
    };

    const_raw_vector<T> getvector(int iStart, int iEnd) const
    {
        #ifndef NO_AP_ASSERT
        if( iEnd>=iStart )
            ap_error::make_assertion(m_iLow<=iStart && iEnd<=m_iHigh);
        #endif
        if( iStart>iEnd )
            return const_raw_vector<T>(0, 0, 1);
        else
            return const_raw_vector<T>(m_Vec+iStart-m_iLow, iEnd-iStart+1, 1);
    };
private:
    T         *m_Vec;
    long      m_iVecSize;
    long      m_iLow, m_iHigh;
};



/********************************************************************
КЛАСС-ШАБЛОН ДВУХМЕРНОГО ДИНАМИЧЕСКОГО МАССИВА


ОПИСАНИЕ ЧЛЕНОВ КЛАССА:

---------------------------------------------------------------------
template_2d_array()

Создание пустого массива.

---------------------------------------------------------------------
~template_2d_array()

Удаление массива. При этом освобождается выделенная под массив память

---------------------------------------------------------------------
template_2d_array(const template_2d_array &rhs)

Создание копии массива. При этом выделяется отдельная область памяти,
в которую копируется содержимое массива-источника

---------------------------------------------------------------------
const template_2d_array& operator=(const template_2d_array &rhs)

Присваивание массива. При этом содержимое массива-приемника удаляется
и освобождается  выделенная под него память,  затем заново выделяется 
отдельная область памяти, в которую копируется содержимое источника.

---------------------------------------------------------------------
      T& operator()(int i1, int i2)
const T& operator()(int i1, int i2) const

Обращение к элементу массива с индексом [i1,i2]

---------------------------------------------------------------------
void setbounds( int iLow1, int iHigh1, int iLow2, int iHigh2 )

Выделение  памяти   под   массив.    При   этом   старое   содержимое
массива   удаляется   и  освобождается  выделенная  под  него  память, 
затем   заново   выделяется   отдельная   область   памяти   размером 
(iHigh1-iLow1+1)*(iHigh2-iLow2+1) элементов.

Нумерация  элементов в новом массиве по первой размерности начинается 
с iLow1 и заканчивается iHigh1, аналогично для второй размерности.

Содержимое нового массива не определено.

---------------------------------------------------------------------
void setcontent( int iLow1, int iHigh1, int iLow2, int iHigh2,
    const T *pContent )

Метод  аналогичен  методу  setbounds()  за тем исключением, что после 
выделения памяти в неё копируется содержимое массива pContent[].

Массив pContent содержит двухмерный массив, записанный построчно, т.е.
первым идет элемент [iLow1, iLow2], затем [iLow1, iLow2+1] и т.д.

---------------------------------------------------------------------
      T* getcontent()
const T* getcontent() const

Метод позволяет получить указатель на содержимое массива. Данные,  на
которые указывает возвращенный указатель, можно изменять, и при  этом
изменится содержимое массива.

---------------------------------------------------------------------
int getlowbound(int iBoundNum)
int gethighbound(int iBoundNum)

Методы используются для  получения  информации  о  нижней  и  верхней
границах массива по размерности с переданным номером.

---------------------------------------------------------------------
raw_vector<T> getcolumn(int iColumn, int iRowStart, int iRowEnd)
const_raw_vector<T> getcolumn(int iColumn, int iRowStart, int iRowEnd) const

Методы используются базовыми подпрограммами линейной алгебры для
получения доступа к внутренней памяти массива. Методы возвращают
объект, содержащий в себе указатель на часть столбца iColumn (начиная
со строки iRowStart и заканчивая строкой iRowEnd).

Параметр iColumn должен быть допустимым номером столбца (т.е. находиться
в пределах выделенной под массив памяти). Если iRowEnd<iRowStart, то
считается, что задан пустой столбец.

---------------------------------------------------------------------
raw_vector<T> getrow(int iRow, int iColumnStart, int iColumnEnd)
const_raw_vector<T> getrow(int iRow, int iColumnStart, int iColumnEnd) const

Методы используются базовыми подпрограммами линейной алгебры для
получения доступа к внутренней памяти массива. Методы возвращают
объект, содержащий в себе указатель на часть строки iRow (начиная
со столбца iColumnStart и заканчивая столбцом iColumnEnd).

Параметр iRow должен быть допустимым номером строки (т.е. находиться
в пределах выделенной под массив памяти). Если iColumnEnd<iColumnStart,
то считается, что задана пустая строка.
********************************************************************/
template<class T>
class template_2d_array
{
public:
    template_2d_array()
    {
        m_Vec=0;
        m_iVecSize=0;
    };

    ~template_2d_array()
    {
        if(m_Vec)
            delete[] m_Vec;
    };

    template_2d_array(const template_2d_array &rhs)
    {
        m_iVecSize = rhs.m_iVecSize;
        m_iLow1 = rhs.m_iLow1;
        m_iLow2 = rhs.m_iLow2;
        m_iHigh1 = rhs.m_iHigh1;
        m_iHigh2 = rhs.m_iHigh2;
        m_iConstOffset = rhs.m_iConstOffset;
        m_iLinearMember = rhs.m_iLinearMember;
        if(rhs.m_Vec)
        {
            m_Vec = new T[m_iVecSize];
            #ifndef UNSAFE_MEM_COPY
            for(int i=0; i<m_iVecSize; i++)
                m_Vec[i] = rhs.m_Vec[i];
            #else
            memcpy(m_Vec, rhs.m_Vec, m_iVecSize*sizeof(T));
            #endif
        }
        else
            m_Vec=0;
    };
    const template_2d_array& operator=(const template_2d_array &rhs)
    {
        if( this==&rhs )
            return *this;

        m_iLow1 = rhs.m_iLow1;
        m_iLow2 = rhs.m_iLow2;
        m_iHigh1 = rhs.m_iHigh1;
        m_iHigh2 = rhs.m_iHigh2;
        m_iConstOffset = rhs.m_iConstOffset;
        m_iLinearMember = rhs.m_iLinearMember;
        m_iVecSize = rhs.m_iVecSize;
        if(m_Vec)
            delete[] m_Vec;
        if(rhs.m_Vec)
        {
            m_Vec = new T[m_iVecSize];
            #ifndef UNSAFE_MEM_COPY
            for(int i=0; i<m_iVecSize; i++)
                m_Vec[i] = rhs.m_Vec[i];
            #else
            memcpy(m_Vec, rhs.m_Vec, m_iVecSize*sizeof(T));
            #endif
        }
        else
            m_Vec=0;
        return *this;
    };

    const T& operator()(int i1, int i2) const
    {
        #ifndef NO_AP_ASSERT
        ap_error::make_assertion(i1>=m_iLow1 && i1<=m_iHigh1);
        ap_error::make_assertion(i2>=m_iLow2 && i2<=m_iHigh2);
        #endif
        return m_Vec[ m_iConstOffset + i2 +i1*m_iLinearMember];
    };

    T& operator()(int i1, int i2)
    {
        #ifndef NO_AP_ASSERT
        ap_error::make_assertion(i1>=m_iLow1 && i1<=m_iHigh1);
        ap_error::make_assertion(i2>=m_iLow2 && i2<=m_iHigh2);
        #endif
        return m_Vec[ m_iConstOffset + i2 +i1*m_iLinearMember];
    };

    void setbounds( int iLow1, int iHigh1, int iLow2, int iHigh2 )
    {
        if(m_Vec)
            delete[] m_Vec;
        m_iVecSize = (iHigh1-iLow1+1)*(iHigh2-iLow2+1);
        m_Vec = new T[m_iVecSize];
        m_iLow1  = iLow1;
        m_iHigh1 = iHigh1;
        m_iLow2  = iLow2;
        m_iHigh2 = iHigh2;
        m_iConstOffset = -m_iLow2-m_iLow1*(m_iHigh2-m_iLow2+1);
        m_iLinearMember = (m_iHigh2-m_iLow2+1);
    };

    void setcontent( int iLow1, int iHigh1, int iLow2, int iHigh2, const T *pContent )
    {
        setbounds(iLow1, iHigh1, iLow2, iHigh2);
        for(int i=0; i<m_iVecSize; i++)
            m_Vec[i]=pContent[i];
    };

    T* getcontent()
    {
        return m_Vec;
    };

    const T* getcontent() const
    {
        return m_Vec;
    };

    int getlowbound(int iBoundNum) const
    {
        return iBoundNum==1 ? m_iLow1 : m_iLow2;
    };

    int gethighbound(int iBoundNum) const
    {
        return iBoundNum==1 ? m_iHigh1 : m_iHigh2;
    };

    raw_vector<T> getcolumn(int iColumn, int iRowStart, int iRowEnd)
    {
        #ifndef NO_AP_ASSERT
        ap_error::make_assertion(m_iLow2<=iColumn && iColumn<=m_iHigh2);
        if( iRowEnd>=iRowStart )
            ap_error::make_assertion(m_iLow1<=iRowStart && iRowEnd<=m_iHigh1);
        #endif
        if( iRowStart>iRowEnd )
            return raw_vector<T>(0, 0, 1);
        else
            return raw_vector<T>(&((*this)(iRowStart, iColumn)), iRowEnd-iRowStart+1, m_iLinearMember);
    };

    raw_vector<T> getrow(int iRow, int iColumnStart, int iColumnEnd)
    {
        #ifndef NO_AP_ASSERT
        ap_error::make_assertion(m_iLow1<=iRow && iRow<=m_iHigh1);
        if( iColumnEnd>=iColumnStart )
            ap_error::make_assertion(m_iLow2<=iColumnStart && iColumnEnd<=m_iHigh2);
        #endif
        if( iColumnStart>iColumnEnd )
            return raw_vector<T>(0, 0, 1);
        else
            return raw_vector<T>(&((*this)(iRow, iColumnStart)), iColumnEnd-iColumnStart+1, 1);
    };

    const_raw_vector<T> getcolumn(int iColumn, int iRowStart, int iRowEnd) const
    {
        #ifndef NO_AP_ASSERT
        ap_error::make_assertion(m_iLow2<=iColumn && iColumn<=m_iHigh2);
        if( iRowEnd>=iRowStart )
            ap_error::make_assertion(m_iLow1<=iRowStart && iRowEnd<=m_iHigh1);
        #endif
        if( iRowStart>iRowEnd )
            return const_raw_vector<T>(0, 0, 1);
        else
            return const_raw_vector<T>(&((*this)(iRowStart, iColumn)), iRowEnd-iRowStart+1, m_iLinearMember);
    };

    const_raw_vector<T> getrow(int iRow, int iColumnStart, int iColumnEnd) const
    {
        #ifndef NO_AP_ASSERT
        ap_error::make_assertion(m_iLow1<=iRow && iRow<=m_iHigh1);
        if( iColumnEnd>=iColumnStart )
            ap_error::make_assertion(m_iLow2<=iColumnStart && iColumnEnd<=m_iHigh2);
        #endif
        if( iColumnStart>iColumnEnd )
            return const_raw_vector<T>(0, 0, 1);
        else
            return const_raw_vector<T>(&((*this)(iRow, iColumnStart)), iColumnEnd-iColumnStart+1, 1);
    };
private:
    T           *m_Vec;
    long        m_iVecSize;
    long        m_iLow1, m_iLow2, m_iHigh1, m_iHigh2;
    long        m_iConstOffset, m_iLinearMember;
};


typedef template_1d_array<int>    integer_1d_array;
typedef template_1d_array<double> real_1d_array;
typedef template_1d_array<int>   boolean_1d_array;
typedef template_2d_array<int>    integer_2d_array;
typedef template_2d_array<double> real_2d_array;
typedef template_2d_array<int>   boolean_2d_array;


/********************************************************************
КОНСТАНТЫ И ФУНКЦИИ, СОВМЕСТИМЫЕ С ALGOPASCAL
********************************************************************/
static double machineepsilon = 5E-16;
static double maxrealnumber = 1E300;
static double minrealnumber = 1E-300;

static double sign(double x)
{
    if( x>0 ) return  1.0;
    if( x<0 ) return -1.0;
    return 0;
}

static double randomreal()
{
    int i = rand();
    while(i==RAND_MAX)
        i =rand();
    return double(i)/double(RAND_MAX);
}

static int randominteger(int maxv)
{  return rand()%maxv; }

static double round(double x)
{ return floor(x+0.5); }

static double trunc(double x)
{ return x>0 ? floor(x) : ceil(x); }

static double pi()
{ return 3.14159265358979323846; }

static double sqr(double x)
{ return x*x; }

static int maxint(int m1, int m2)
{
    return m1>m2 ? m1 : m2;
}

static int minint(int m1, int m2)
{
    return m1>m2 ? m2 : m1;
}

static double maxreal(double m1, double m2)
{
    return m1>m2 ? m1 : m2;
}

static double minreal(double m1, double m2)
{
    return m1>m2 ? m2 : m1;
}

};//namespace ap


#endif
