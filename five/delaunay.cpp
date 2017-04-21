const double EP = 0.00000001;

// 点结构
typedef struct PT{
    double x;
    double y;
}PT;

//线结构
typedef struct SEGMENT
{
    PT* ptStart;
    PT* ptEnd;
}SEGMENT;

//三角形类
class TRIANGLE
{
public:
    TRIANGLE(PT* pt1, PT* pt2, PT* pt3)
    {
        //确保三个点x坐标升序排序
        PT* temp;
        if(pt1->x > pt2->x){temp = pt1; pt1 = pt2; pt2 = temp;}
        if(pt3->x < pt1->x  && pt3->x )
        {
            temp = pt3; pt3 = pt1; pt1 = temp;
            if(pt3->x < pt2->x){temp=pt3;pt3=pt2;pt2=temp;}
        }
        if(pt3->x < pt2->x){temp=pt3;pt3=pt2;pt2=temp;}

        _ptFirst = pt1;
        _ptSecond = pt2;
        _ptThird = pt3;

        InitData();

    }
    PT* GetFirstPt(){return _ptFirst;}
    PT* GetSecondPt(){return _ptSecond;}
    PT* GetThirdPt(){return _ptThird;}

    SEGMENT GetSegment1(){SEGMENT seg = {_ptFirst,_ptSecond};return seg;}
    SEGMENT GetSegment2(){SEGMENT seg = {_ptFirst,_ptThird};return seg;}
    SEGMENT GetSegment3(){SEGMENT seg = {_ptSecond,_ptThird};return seg;}

    BOOL IsTriangleSame(TRIANGLE& tri)
    {
        if(tri.GetFirstPt() == _ptFirst &&
            tri.GetSecondPt() == _ptSecond &&
            tri.GetThirdPt() == _ptThird)
            return TRUE;

        return FALSE;
    }

    BOOL IsPtInCircle(PT* pt)
    {
        double offsetx = pt->x - _ptCenter.x;
        double offsety = pt->y - _ptCenter.y;
        if(sqrt(offsetx*offsetx + offsety*offsety) <= _Radius)
            return TRUE;
    
        return FALSE;
                        
    }
protected:

  // 获取三角形外接圆中心点及半径
    void InitData()
    {
        double x0 = _ptFirst->x;
        double y0 = _ptFirst->y;
        
        double x1 = _ptSecond->x;
        double y1 = _ptSecond->y;
        
        double x2 = _ptThird->x;
        double y2 = _ptThird->y;
        
        
        double y10 = y1 - y0;
        double y21 = y2 - y1;
        
        bool b21zero = y21 > -REAL_EPSILON && y21 < REAL_EPSILON;
        
        if (y10 > -REAL_EPSILON && y10 < REAL_EPSILON)
        {
            if (b21zero)  
            {
                if (x1 > x0)
                {
                    if (x2 > x1) x1 = x2;
                }
                else
                {
                    if (x2 < x0) x0 = x2;
                }
                _ptCenter.x = (x0 + x1) * .5F;
                _ptCenter.y = y0;
            }
            else  
            {
                double m1 = - (x2 - x1) / y21;
                
                double mx1 = (x1 + x2) * .5F;
                double my1 = (y1 + y2) * .5F;
                
                _ptCenter.x = (x0 + x1) * .5F;
                _ptCenter.y = m1 * (_ptCenter.x - mx1) + my1;
            }
        }
        else if (b21zero)  
        {
            double m0 = - (x1 - x0) / y10;
            
            double mx0 = (x0 + x1) * .5F;
            double my0 = (y0 + y1) * .5F;
            
            _ptCenter.x = (x1 + x2) * .5F;
            _ptCenter.y = m0 * (_ptCenter.x - mx0) + my0;
        }
        else  
        {
            double m0 = - (x1 - x0) / y10;
            double m1 = - (x2 - x1) / y21;
            
            double mx0 = (x0 + x1) * .5F;
            double my0 = (y0 + y1) * .5F;
            
            double mx1 = (x1 + x2) * .5F;
            double my1 = (y1 + y2) * .5F;
            
            _ptCenter.x = (m0 * mx0 - m1 * mx1 + my1 - my0) / (m0 - m1);
            _ptCenter.y = m0 * (_ptCenter.x - mx0) + my0;
        }
        
        double dx = x0 - _ptCenter.x;
        double dy = y0 - _ptCenter.y;
        
        _Radius2 = dx * dx + dy * dy;    // the radius of the circumcircle, squared
        _Radius = (double) sqrt(_Radius2);    // the proper radius
      
        _Radius2 *= 1.000001f;
}
private:
    PT* _ptFirst;
    PT* _ptSecond;
    PT* _ptThird;

    PT            _ptCenter;
    double        _Radius2;
    double        _Radius;
};

class CMyDelaunay  
{
public:

// 通过样本点集，生成三角形集，由vecTriangleWork输出
    void BuildDelaunayEx(vector<PT>& vecPT, vector<TRIANGLE>& vecTriangleWork);

    BOOL IsPtsBuildTriangle(PT* pt1, PT* pt2, PT* pt3);
public:
    CMyDelaunay();
    virtual ~CMyDelaunay();

};

//判断三个点能否组成一个三角形
BOOL CMyDelaunay::IsPtsBuildTriangle(PT* pt1, PT* pt2, PT* pt3)
{
    double offset_x1 = pt2->x - pt1->x;
    double offset_x2 = pt3->x - pt2->x;
    double offset_y1 = pt2->y - pt1->y;
    double offset_y2 = pt3->y - pt2->y;

    if((fabs(offset_x1) < EP) && (fabs(offset_x2) < EP))    //竖直
    {
        return FALSE;
    }
    
    if(fabs(offset_x1) > EP && fabs(offset_x2) > EP)
    {
        if(fabs(offset_y1/offset_x1 - offset_y2/offset_x2) < EP)
            return FALSE;
    }
    
    return TRUE;
}

void CMyDelaunay::BuildDelaunayEx(vector<PT>& vecPT, vector<TRIANGLE>& vecTriangleWork)
{
    int nSize = vecPT.size();
    if(nSize < 3)
        return;
    
    for(int i = 0; i < nSize - 2; ++i)
    {
        for(int j = i + 1;  j < nSize - 1; ++j)
        {
            for(int k = j + 1; k < nSize; ++k)
            {
                PT* pt1 = &vecPT[i];
                PT* pt2 = &vecPT[j];
                PT* pt3 = &vecPT[k];
                BOOL bFind = TRUE;
                for(int m = 0; m < nSize; ++m)
                {
                    PT* pt = &vecPT[m];
                    if(pt != pt1 && pt != pt2 && pt != pt3 && IsPtsBuildTriangle(pt1, pt2, pt3))
                    {
                        TRIANGLE tri(pt1, pt2, pt3);
                        if(tri.IsPtInCircle(pt))
                        {
                            bFind = FALSE;
                            break;
                        }
                    }
                }
                if(bFind)
                {  
                    vecTriangleWork.push_back(tri);  
                }

            }
        }
    }
}