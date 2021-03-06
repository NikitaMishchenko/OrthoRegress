#ifndef ORTHOGONAL_POLYNOMIAL_CURVE_FITTING_H_INCLUDED
#define ORTHOGONAL_POLYNOMIAL_CURVE_FITTING_H_INCLUDED

///reference https://ru.stackoverflow.com/questions/383716/%D0%9C%D0%9D%D0%9A-%D0%B0%D0%BB%D0%B3%D0%BE%D1%80%D0%B8%D1%82%D0%BC

///DOC here: https://www.google.com.ua/url?sa=t&rct=j&q=&esrc=s&source=web&cd=8&cad=rja&uact=8&ved=0ahUKEwjYw4eXg9HSAhXPKiwKHRJBBU0QFghXMAc&url=http%3A%2F%2Fjeffareid.net%2Fmisc%2Fopls.rtf&usg=AFQjCNHuzP-Q0c5BABak1M7NhxKN29Gr0w&sig2=N6J8HAijGQOuA6qJrOJc3A&bvm=bv.149397726,d.bGg
#include <iostream>
#include <iomanip>
#include <cmath>

#include <vector>

class Ortho_GLS
{
protected:

    std::vector<std::vector<double>> P;
public:
    enum MODES
    {
        CONSOLE_REPORT_MAXINFO,
        CONSOLE_REPORT_EACH_STEP,
        CONSOLE_SILENCE
    };
    MODES Regime = MODES::CONSOLE_SILENCE;


    void f_QjSj(int j, const std::vector<double> &x, std::vector<double> &Qj, std::vector<double> &Sj)
    {//std::cout << "void f_QjSj()\n";
        if(j == 0){
            Qj[j] = Sj[j] = 0.0;
            for(int i = 0; i < x.size(); ++i){
                Qj[j] += x[i];
                Sj[j] += pow(P[j][i],2);///Sj = m+1 always//if j == 0
            }//std::cout << "f_QjSj j = " << j << std::endl;
        }

        if(j == 1){
            Qj[j] = Sj[j] = 0.0;
            for(int i = 0; i < x.size(); ++i){//std::cout << i << std::endl;
                Qj[j] += x[i]*pow(P[j][i],2);
                Sj[j] += pow(P[j][i],2);
            }//std::cout << "f_QjSj j = " << j << std::endl;
        }

        if(j > 1){
            Qj[j] = Sj[j] = 0.0;
            for(int i = 0; i < x.size(); ++i){
                Qj[j] += x[i]*pow(P[j][i],2);
                Sj[j] += pow(P[j][i],2);
            }//std::cout << "f_QjSj j = " << j << std::endl;
        }
    }

    void f_QjSj(int j, const std::vector<double> &x, const std::vector<double> &w, std::vector<double> &Qj, std::vector<double> &Sj)
    {//std::cout << "void f_QjSj()\n";
        if(j == 0){
            Qj[j] = Sj[j] = 0.0;
            for(int i = 0; i < x.size(); ++i){
                Qj[j] += w[i]*x[i];
                Sj[j] += w[i]*pow(P[j][i],2);///Sj = m+1 always//if j == 0
            }//std::cout << "f_QjSj j = " << j << std::endl;
        }

        if(j == 1){
            Qj[j] = Sj[j] = 0.0;
            for(int i = 0; i < x.size(); ++i){//std::cout << i << std::endl;
                Qj[j] += w[i]*x[i]*pow(P[j][i],2);
                Sj[j] += w[i]*pow(P[j][i],2);
            }//std::cout << "f_QjSj j = " << j << std::endl;
        }

        if(j > 1){
            Qj[j] = Sj[j] = 0.0;
            for(int i = 0; i < x.size(); ++i){
                Qj[j] += w[i]*x[i]*pow(P[j][i],2);
                Sj[j] += w[i]*pow(P[j][i],2);
            }//std::cout << "f_QjSj j = " << j << std::endl;
        }
    }

    void orthogonal_polynom_P(int m, const std::vector<double> &x ,std::vector<double> &_S, std::vector<double> &_Q)
    {
        int j;
        j = 0;
        if(m >= j){
            for(int i = 0; i < x.size(); ++i){
                P[0][i] = 1.0;
            this->f_QjSj(j, x, _Q, _S);//std::cout << "Q = " << _Q[j] << "\tS = " << _S[j] << std::endl;
            }
        }

        j = 1;
        if(m >= j){
            for(int i = 0; i < x.size(); ++i){
                P[j][i] = x[i] - _Q[j-1]/_S[j-1];
            this->f_QjSj(j, x, _Q, _S);//std::cout << "Q = " << _Q[j] << "\tS = " << _S[j] << std::endl;
            }
        }

        j=2;
        if(j >= 2 && j <= m+1){
            for(; j <= m+1; j++){
                for(int i = 0; i < x.size(); ++i){
                    P[j][i] = (x[i] - _Q[j-1]/_S[j-1])*P[j-1][i] - _S[j-1]/_S[j-2] * P[j-2][i];
                this->f_QjSj(j, x, _Q, _S);//std::cout << "Q = " << _Q[j] << "\tS = " << _S[j] << std::endl;
                }
            }
        }
    }

    void orthogonal_polynom_P(int m, const std::vector<double> &x, const std::vector<double> &w, std::vector<double> &_S, std::vector<double> &_Q)
    {
        int j;
        j = 0;
        if(m >= j){
            for(int i = 0; i < x.size(); ++i){
                P[0][i] = 1.0;
            this->f_QjSj(j, x, _Q, _S);//std::cout << "Q = " << _Q[j] << "\tS = " << _S[j] << std::endl;
            }
        }

        j = 1;
        if(m >= j){
            for(int i = 0; i < x.size(); ++i){
                P[j][i] = x[i] - _Q[j-1]/_S[j-1];
            this->f_QjSj(j, x, _Q, _S);//std::cout << "Q = " << _Q[j] << "\tS = " << _S[j] << std::endl;
            }
        }

        j=2;
        if(j >= 2 && j <= m+1){
            for(; j <= m+1; j++){
                for(int i = 0; i < x.size(); ++i){
                    P[j][i] = (x[i] - _Q[j-1]/_S[j-1])*P[j-1][i] - _S[j-1]/_S[j-2] * P[j-2][i];
                this->f_QjSj(j, x, _Q, _S);//std::cout << "Q = " << _Q[j] << "\tS = " << _S[j] << std::endl;
                }
            }
        }
    }

    void recurrent_P(int m, const std::vector<double> &x, const std::vector<double> &y, std::vector<double> &polynom)
    {
    if(Regime != MODES::CONSOLE_SILENCE)
    {
        std::cout << "void recurrent_P()\n";

        std::cout << "Array for polynom\n";
        for(auto _x : x)
            std::cout << _x << "\  ";
                std::cout << std::endl;
        for(auto _y : y)
            std::cout << _y << "\  ";
                std::cout << std::endl;
        std::cout << "Degree of polynom = " << m << std::endl;
    }

        std::vector<double> _S(m+1);
        std::vector<double> _Q(m+1);
        std::vector<std::vector<double>> nP(x.size(), std::vector<double> (m+1));
        P = nP;

        orthogonal_polynom_P(m, x, _S, _Q);
    if(Regime != MODES::CONSOLE_SILENCE){
        std::cout << "SUM(P^2) S:\n";
        for(int i = 0; i < _S.size(); i++)
            std::cout << _S[i] << "\t";
        std::cout << std::endl;
        std::cout << "SUM(xP^2) Q:\n";
        for(int i = 0; i < _Q.size(); i++)
            std::cout << _Q[i] << "\t";
        std::cout << std::endl << std::endl;

        std::cout << "POLYNOM COEFFICIENTS P:\n";
        for(int i = 0; i < m+1; i++)
        {
            std::cout << "i->" << i << "\t";
            for(int j = 0; j < x.size(); j++)
                std::cout << P[i][j] << "\t";
            std::cout << std::endl;
        }
        std::cout << std::endl << std::endl;
    }

        std::vector<std::vector<double>> c(m+2, std::vector<double> (m+2));
        orthogonal_polynom_coefficients(m, _S, _Q, c);

    if(Regime != MODES::CONSOLE_SILENCE)
    {
        std::cout << "POLYNOM COEFF c:\n";
        for(int i = 0; i < c.size(); i++)
        {
            std::cout << "i->" << i << "\t";
            for(int j = 0; j < c[0].size(); j++)
                std::cout << c[i][j] << "\t";
            std::cout << std::endl;
        }
        std::cout << std::endl << std::endl;
    }
        std::vector<double> a(m+1);
        std::vector<double> b(m+1);
        polynom_coefficeints(m, y, _S, c, b, a);
            polynom = a;

    if(Regime != MODES::CONSOLE_SILENCE)
    {
        std::cout << "B:\n";
        for(int i = 0; i < b.size(); i++)
            std::cout << std::setprecision(12) << b[i] << "\t";
        std::cout << std::endl;

        std::cout << "A:\n";
        for(int i = 0; i < a.size(); i++)
            std::cout << std::setprecision(12) << a[i] << "\t";
        std::cout << std::endl;
    }
    };

    void orthogonal_polynom_coefficients(int m, const std::vector<double> S, const std::vector<double> Q, std::vector<std::vector<double>> &c)
    {
        ///k height, j width
        c[0][0] = 1.0;///k = 0
        c[1][0] = -Q[0]/S[0];   c[1][1] = 1.0; ///k = 1

        for(int k = 2; k < c.size(); k++)
            for(int j = 0; j <= k; ++j)
                if(k != j)
                    c[k][j] = (j == 0 ? 0.0 : c[k-1][j-1]) - Q[k-1]/S[k-1]*c[k-1][j] - (j == k ? 0.0 : S[k-1]/S[k-2]*c[k-2][j]);
                else
                    c[k][j] = 1.0;
    }

    void polynom_coefficeints(int m, const std::vector<double> y, const std::vector<double> S, const std::vector<std::vector<double>> &c,
                              std::vector<double> &b, std::vector<double> &a)///b.size == m
    {
        ///g(x) = SUM(b0*p0(x) + b1*p1(x) + ... + bm*pm(x))
        for(int i = 0; i < b.size(); ++i)
            b[i] = 0;

        for(int k = 0; k < b.size(); ++k)
            for(int i = 0; i < y.size(); ++i)
                b[k] += y[i]*P[k][i]/S[k];

        if(Regime != MODES::CONSOLE_SILENCE)
        {
            std::cout << "b[] calculated\n";
        }

        ///g(x) = SUM(a0 + a1*x + a2*x^2 + ... + am*x^m)
        for(int i = 0; i < a.size(); ++i)
            a[i] = 0.0;

        for(int j = 0; j < c.size(); ++j)
            for(int k = j; k < c[0].size()-1; ++k)
                    a[j] += b[k]*c[k][j];

        if(Regime != MODES::CONSOLE_SILENCE)
        {
            std::cout << "a[] calculated\n\n";
        }
    }

    void Test()
    {std::cout << "Algorithm TEST\n";
        Regime = MODES::CONSOLE_REPORT_MAXINFO;

        ///input
        std::vector<double> x = {0, 3.36588, 3.63719, 0.56448, -3.02721, -3.8357, -1.11766, 2.62795, 3.95743, 1.64847};
        std::vector<double> y = {3.9561, 74.84479, 89.44289, 6.46668, -14.53888, -34.55881, 1.70531, 43.80101, 109.1294, 18.81613};
        int m = 3;

        std::vector<double> _S(m+1);
        std::vector<double> _Q(m+1);
        std::vector<std::vector<double>> nP(x.size(), std::vector<double> (m+1));
            P = nP;

        std::cout << "orthogonal_polynom_P() started\n";
        orthogonal_polynom_P(m, x, _S, _Q);
            std::cout << "_S:";
            std::cout << _S[0] << "\t" << fabs(_S[0] - 10) << "\n";
            std::cout << _S[1] << "\t" << fabs(_S[1] - 69.17098425001) << "\n";
            std::cout << _S[2] << "\t" << fabs(_S[2] - 328.7768235086) << "\n";
            std::cout << _S[3] << "\t" << fabs(_S[3] - 962.75401849569) << "\n";
                std::cout << std::endl;
            std::cout << "_Q:";
            std::cout << _Q[0] << "\t" << fabs(_Q[0] - 7.82083) << "\n";
            std::cout << _Q[1] << "\t" << fabs(_Q[1] - (-27.512890311149)) << "\n";
            std::cout << _Q[2] << "\t" << fabs(_Q[2] - (-1.7129366948718)) << "\n";
            std::cout << _Q[3] << "\t" << fabs(_Q[3] - 250.39738237352) << "\n";
                        std::cout << std::endl;

        std::cout << "orthogonal_polynom_coefficients() started\n";
        std::vector<std::vector<double>> c(m+2, std::vector<double> (m+2));
        orthogonal_polynom_coefficients(m, _S, _Q, c);
        //0 => [1]
        //1 => [-0.782083, 1]
        //2 => [-7.2281734230875, -0.38433110143359, 1]
        //3 => [3.6796621840027, -11.983278954678, -0.37912107270775, 1]
        //4 => [20.209167857339, 7.9217601883617, -14.812965853531, -0.63920555697172, 1]
        std::cout << "err count started\n";
        std::vector<double> err_orthogonal_polynom;
        err_orthogonal_polynom.push_back(abs(1.0-c[0][0]));

        err_orthogonal_polynom.push_back(abs(-0.782083-c[1][0]));
        err_orthogonal_polynom.push_back(abs(1.0-c[1][1]));

        err_orthogonal_polynom.push_back(abs(-7.2281734230875-c[2][0]));
        err_orthogonal_polynom.push_back(abs(-0.38433110143359-c[2][1]));
        err_orthogonal_polynom.push_back(abs(1.0-c[2][2]));

        err_orthogonal_polynom.push_back(abs(3.6796621840027-c[3][0]));
        err_orthogonal_polynom.push_back(abs(-11.983278954678-c[3][1]));
        err_orthogonal_polynom.push_back(abs(-0.37912107270775-c[3][2]));
        err_orthogonal_polynom.push_back(abs(1.0-c[3][3]));

        err_orthogonal_polynom.push_back(abs(20.209167857339-c[4][0]));
        err_orthogonal_polynom.push_back(abs(7.9217601883617-c[4][1]));
        err_orthogonal_polynom.push_back(abs(-14.812965853531-c[4][2]));
        err_orthogonal_polynom.push_back(abs(-0.63920555697172-c[4][3]));
        err_orthogonal_polynom.push_back(abs(1.0-c[4][4]));


        std::cout << "ERRs orthogonal_polynom_coefficients(m, _S, _Q, c)\n";
        for(auto err : err_orthogonal_polynom)
            std::cout << std::setprecision(9) << err << "\t";
                std::cout << std::endl << std::endl;


        std::vector<double> a(m+1);
        std::vector<double> b(m+1);
        polynom_coefficeints(m, y, _S, c, b, a);

        //b => [29.906462, 15.897655646369, 2.3791287087584, 1.000001267701]
        //a => [3.9560877250835, 2.9999883433859, 2.0000071554385, 1.000001267701]
        std::cout << "ERRs polynom_coefficeints(m, y, _S, c, b, a)\n";
        std::cout << "b\n";
        std::vector<double> err_b;
            err_b.push_back(abs(29.906462-b[0]));
            err_b.push_back(abs(15.897655646369-b[1]));
            err_b.push_back(abs(2.3791287087584-b[2]));
            err_b.push_back(abs(1.000001267701-b[3]));
        for(auto err : err_b)
            std::cout << std::setprecision(9) << err << "\t";
                std::cout << std::endl << std::endl;

        std::cout << "a\n";
        std::vector<double> err_a;
            err_a.push_back(abs(3.9560877250835-a[0]));
            err_a.push_back(abs(2.9999883433859-a[1]));
            err_a.push_back(abs(2.0000071554385-a[2]));
            err_a.push_back(abs(1.000001267701-a[3]));
        for(auto err : err_a)
            std::cout << std::setprecision(9) << err << "\t";
                std::cout << std::endl << std::endl;
    }
};



#endif // ORTHOGONAL_POLYNOMIAL_CURVE_FITTING_H_INCLUDED
