#ifndef LOTKA_VOLTERRA_H
#define LOTKA_VOLTERRA_H

#include <fstream>
#include <vector>
#include <iostream>

class Lotka_Volterra {

private:

	double alpha;
	double beta;
	double gamma;
	double delta;

	double x0;
	double y0;

    std::pair<double, double> solveNewton(double x_prev, double y_prev, double h) {

        double x_k = x_prev;
        double y_k = y_prev;

        //x_(k+1) = x_(k) - [DF(x_k)^(-1)]F(x_k)

        //(x_(k+1), y_(k+1)) = (x_k, y_k) - [DF^-1] * (f1(x, y), f2(x, y))

        //(x_(k+1), y_(k+1)) - (x_k, y_k) = -[DF(x_k)^(-1)] * (f1, f2)
        
        //(dx, dy) = -[DF(x_k)^(-1)] * (f1, f2)

        //[DF(x_k)](dx, dy)= -(f1, f2)

        //Ax = f, Крамер

        for (int i = 0; i < 100; ++i) {

            double f1 = x_k - x_prev - 0.5 * h * (alpha * x_k - beta * x_k * y_k + alpha * x_prev - beta * x_prev * y_prev); 
            double f2 = y_k - y_prev - 0.5 * h * (delta * x_k * y_k - gamma * y_k + delta * x_prev * y_prev - gamma * y_prev);

            double df1_dx = 1 - 0.5 * h * (alpha - beta * y_k);
            double df1_dy = 0.5 * h * beta * x_k;
            double df2_dx = 0.5 * h * delta * y_k;
            double df2_dy = 1 - 0.5 * h * (-gamma + delta * x_k);

            double det = df1_dx * df2_dy - df1_dy * df2_dx;



            double dx = (f1 * df2_dy - f2 * df1_dy) / det;
            double dy = (df1_dx * f2 - df2_dx * f1) / det;

            x_k -= dx;
            y_k -= dy;

            if (std::abs(dx) < 1e-6 && std::abs(dy) < 1e-6) {
                 break;
            }

        }
        return std::make_pair(x_k, y_k);

    }

    std::pair<double, double> dt(double x, double y) {

        double dxdt = alpha * x - beta * x * y;
        double dydt = delta * x * y - gamma * y;

        return std::make_pair(dxdt, dydt);

    }

public:

	Lotka_Volterra(double a, double b, double g, double d, double x, double y) : alpha{ a }, beta{ b }, gamma{ g }, delta{ d }, x0{ x }, y0{ y } {};


    void trapezoidMethod(double h, double t_end) {

        std::vector<std::pair<double, double>> results;

        double t = 0.0;

        double x = x0;
        double y = y0;

        results.push_back({ t, x });
        results.push_back({ t, y });

        while (t < t_end) {

            std::pair<double, double> next_point = solveNewton(x, y, h);

            t += h;

            x = next_point.first;
            y = next_point.second;

            results.push_back(std::make_pair(t, x));
            results.push_back(std::make_pair(t, y));

        }

        std::ofstream ofs("xy_slv.txt");

        for (int i = 0; i < results.size(); i += 2) {
            ofs << results[i].second << " "
                << results[i + 1].second << '\n';
        }

        ofs.close(); ofs.clear();
        ofs.open("ty_slv.txt");
        for (int i = 1; i < results.size(); i += 2) {
            ofs << results[i].first << " " << results[i].second << '\n';
        }
        ofs.close(); ofs.clear();
        ofs.open("tx_slv.txt");
        for (int i = 0; i < results.size(); i += 2) {
            ofs << results[i].first << " " << results[i].second << '\n';
        }
        ofs.close(); ofs.clear();
    }

    void rungeKutta4(double h, double t_end) {

        double t = 0.0;

        double x = x0;
        double y = y0;

        std::vector<std::pair<double, double>> results;

        while (t < t_end) {

            results.emplace_back(t, x);
            results.emplace_back(t, y);

            std::pair<double, double> k1 = dt(x, y);
            std::pair<double, double> k2 = dt(x + h * k1.first / 2, y + h * k1.second / 2);
            std::pair<double, double> k3 = dt(x + h * k2.first / 2, y + h * k2.second / 2);
            std::pair<double, double> k4 = dt(x + h * k3.first, y + h * k3.second);

            x += h / 6 * (k1.first + 2 * k2.first + 2 * k3.first + k4.first);
            y += h / 6 * (k1.second + 2 * k2.second + 2 * k3.second + k4.second);

            t += h;

        }


        std::ofstream ofs("xy_rk4.txt");
        for (int i = 0; i < results.size(); i += 2) {
            ofs << results[i].second << " " 
                << results[i + 1].second << '\n';
        }
        ofs.close(); ofs.clear();
        ofs.open("ty_rk4.txt");
        for (int i = 1; i < results.size(); i += 2) {
            ofs << results[i].first << " " << results[i].second << '\n';
        }
        ofs.close(); ofs.clear();
        ofs.open("tx_rk4.txt");
        for (int i = 0; i < results.size(); i += 2) {
            ofs << results[i].first << " " << results[i].second << '\n';
        }
        ofs.close(); ofs.clear();
    }

    void eulerMethod(double h, double t_end) {

        double t = 0.0;

        double x = x0;
        double y = y0;

        std::vector<std::pair<double, double>> results;

        while (t < t_end) {

            results.push_back(std::make_pair(t, x));
            results.push_back(std::make_pair(t, y));

            std::pair<double, double> d = dt(x, y);

            x += h * d.first;
            y += h * d.second;

            t += h;

        }

        std::ofstream ofs("xy.txt");

        for (int i = 0; i < results.size(); i += 2) {
            ofs << results[i].second << " "
                << results[i + 1].second << '\n';
        }
        ofs.close(); ofs.clear();
        ofs.open("ty.txt");
        for (int i = 1; i < results.size(); i += 2) {
            ofs << results[i].first << " " << results[i].second << '\n';
        }
        ofs.close(); ofs.clear();
        ofs.open("tx.txt");
        for (int i = 0; i < results.size(); i += 2) {
            ofs << results[i].first << " " << results[i].second << '\n';
        }
        ofs.close(); ofs.clear();
    }

	~Lotka_Volterra() {};

};

#endif
#pragma once
