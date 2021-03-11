//
// Created by Report on 2021/2/6.
//
#include <vector>
#include <pybind11/pybind11.h>
#include<pybind11/stl.h>
#include <iostream>
#include <pybind11/eigen.h>

using namespace std;

namespace py = pybind11;
using namespace Eigen;

vector<std::array<float, 3>> cubic() {
    vector<std::array<float, 3>> kpcubic = { /* NOLINT */
            std::array<float, 3>({0.0, 0.0, 0.0}),
            std::array<float, 3>({0.0, 0.5, 0.0}),
            std::array<float, 3>({0.5, 0.5, 0.5}),
            std::array<float, 3>({0.5, 0.5, 0.0}),
    };
    return kpcubic;
}

vector<std::array<float, 3>> fcc() {
    vector<std::array<float, 3>> kpfcc = { /* NOLINT */
            std::array<float, 3>({0.0, 0.0, 0.0}),
            std::array<float, 3>({3.0 / 8.0, 3.0 / 8.0, 3.0 / 4.0}),
            std::array<float, 3>({0.5, 0.5, 0.5}),
            std::array<float, 3>({5.0 / 8.0, 1.0 / 4.0, 5.0 / 8.0}),
            std::array<float, 3>({0.5, 1.0 / 4.0, 3.0 / 4.0}),
            std::array<float, 3>({0.5, 0.0, 0.5}),
    };
    return kpfcc;
}

vector<std::array<float, 3>> bcc() {
    vector<std::array<float, 3>> kpbcc = { /* NOLINT */
            std::array<float, 3>({0.0, 0.0, 0.0}),
            std::array<float, 3>({0.5, -0.5, 0.5}),
            std::array<float, 3>({0.25, 0.25, 0.25}),
            std::array<float, 3>({0.0, 0.0, 0.5}),
    };
    return kpbcc;
}

vector<std::array<float, 3>> tet() {
    vector<std::array<float, 3>> kptet = {
            std::array<float, 3>({0.0, 0.0, 0.0}),
            std::array<float, 3>({0.5, 0.5, 0.5}),
            std::array<float, 3>({0.5, 0.5, 0.0}),
            std::array<float, 3>({0.0, 0.5, 0.5}),
            std::array<float, 3>({0.0, 0.5, 0.0}),
            std::array<float, 3>({0.0, 0.0, 0.5}),
    };
    return kptet;
}


vector<std::array<float, 3>> bctet1(float c, float a) {
    auto eta = (float) ((1 + pow(c, 2) / pow(a, 2)) / 4.0);

    vector<std::array<float, 3>> kp = {
            std::array<float, 3>({0.0, 0.0, 0.0}),
            std::array<float, 3>({-0.5, 0.5, 0.5}),
            std::array<float, 3>({0.0, 0.5, 0.0}),
            std::array<float, 3>({0.25, 0.25, 0.25}),
            std::array<float, 3>({0.0, 0.0, 0.5}),
            std::array<float, 3>({eta, eta, -eta}),
            std::array<float, 3>({-eta, 1 - eta, eta}),
    };

    return kp;
}


vector<std::array<float, 3>> bctet2(float c, float a) {
    auto eta = (float) ((1 + pow(a, 2) / pow(c, 2)) / 4.0);
    auto zeta = (float) (pow(a, 2) / (2 * pow(c, 2)));

    vector<std::array<float, 3>> kp = {
            std::array<float, 3>({0.0, 0.0, 0.0}),
            std::array<float, 3>({0.0, 0.5, 0.0}),
            std::array<float, 3>({0.25, 0.25, 0.25}),
            std::array<float, 3>({-eta, eta, eta}),
            std::array<float, 3>({eta, 1 - eta, -eta}),
            std::array<float, 3>({0.0, 0.0, 0.5}),
            std::array<float, 3>({-zeta, zeta, 0.5}),
            std::array<float, 3>({0.5, 0.5, -zeta}),
            std::array<float, 3>({0.5, 0.5, -0.5}),
    };

    return kp;
}


vector<std::array<float, 3>> orc() {

    vector<std::array<float, 3>> kp = {
            std::array<float, 3>({0.0, 0.0, 0.0}),
            std::array<float, 3>({0.5, 0.5, 0.5}),
            std::array<float, 3>({0.5, 0.5, 0.0}),
            std::array<float, 3>({0.0, 0.5, 0.5}),
            std::array<float, 3>({0.5, 0.0, 0.5}),
            std::array<float, 3>({0.5, 0.0, 0.0}),
            std::array<float, 3>({0.0, 0.5, 0.0}),
            std::array<float, 3>({0.0, 0.0, 0.5}),
    };

    return kp;
}

vector<std::array<float, 3>> orcf1(float a, float b, float c) {


    auto zeta = (float) ((1 + pow(a, 2) / pow(b, 2) - pow(a, 2) / pow(c, 2)) / 4);
    auto eta = (float) ((1 + pow(a, 2) / pow(b, 2) + pow(a, 2) / pow(c, 2)) / 4);

    vector<std::array<float, 3>> kp = {
            std::array<float, 3>({0.0, 0.0, 0.0}),
            std::array<float, 3>({0.5, static_cast<float>(0.5 + zeta), zeta}),
            std::array<float, 3>({0.5, static_cast<float>(0.5 - zeta), 1 - zeta}),
            std::array<float, 3>({0.5, 0.5, 0.5}),
            std::array<float, 3>({1, 0.5, 0.5}),
            std::array<float, 3>({0.0, eta, eta}),
            std::array<float, 3>({1, 1 - eta, 1 - eta}),
            std::array<float, 3>({0.5, 0.0, 0.5}),
            std::array<float, 3>({0.5, 0.5, 0.0}),
    };

    return kp;
}


vector<std::array<float, 3>> orcf2(float a, float b, float c) {


    auto phi = (float) ((1 + pow(c, 2) / pow(b, 2) - pow(c, 2) / pow(a, 2)) / 4);
    auto eta = (float) ((1 + pow(a, 2) / pow(b, 2) - pow(a, 2) / pow(c, 2)) / 4);
    auto delta = (float) ((1 + pow(b, 2) / pow(a, 2) - pow(b, 2) / pow(c, 2)) / 4);
    vector<std::array<float, 3>> kp = {
            std::array<float, 3>({0.0, 0.0, 0.0}),
            std::array<float, 3>({0.5, static_cast<float>(0.5 - eta), 1 - eta}),
            std::array<float, 3>({0.5, static_cast<float>(0.5 + eta), eta}),
            std::array<float, 3>({static_cast<float>(0.5 - delta), 0.5, 1 - delta}),
            std::array<float, 3>({static_cast<float>(0.5 + delta), 0.5, delta}),
            std::array<float, 3>({0.5, 0.5, 0.5}),
            std::array<float, 3>({1 - phi, static_cast<float>(0.5 - phi), 0.5}),
            std::array<float, 3>({phi, static_cast<float>(0.5 + phi), 0.5}),
            std::array<float, 3>({0.0, 0.5, 0.5}),
            std::array<float, 3>({0.5, 0.0, 0.5}),
            std::array<float, 3>({0.5, 0.5, 0.0}),
    };

    return kp;
}


vector<std::array<float, 3>> orcf3(float a, float b, float c) {


    auto zeta = float((1 + pow(a, 2) / pow(b, 2) - pow(a, 2) / pow(c, 2)) / 4);
    auto eta = float((1 + pow(a, 2) / pow(b, 2) + pow(a, 2) / pow(c, 2)) / 4);
    vector<std::array<float, 3>> kp = {
            std::array<float, 3>({0.0, 0.0, 0.0}),
            std::array<float, 3>({0.5, static_cast<float>(0.5 + zeta), zeta}),
            std::array<float, 3>({0.5, static_cast<float>(0.5 - zeta), 1 - zeta}),
            std::array<float, 3>({0.5, 0.5, 0.5}),
            std::array<float, 3>({1, 0.5, 0.5}),
            std::array<float, 3>({0.0, eta, eta}),
            std::array<float, 3>({1, 1 - eta, 1 - eta}),
            std::array<float, 3>({0.5, 0.0, 0.5}),
            std::array<float, 3>({0.5, 0.5, 0.0}),
    };

    return kp;
}


vector<std::array<float, 3>> orci(float a, float b, float c) {


    auto zeta = (float) ((1 + pow(a, 2) / pow(c, 2)) / 4);
    auto eta = (float) ((1 + pow(b, 2) / pow(c, 2)) / 4);
    auto delta = (float) ((pow(b, 2) - pow(a, 2)) / (4 * pow(c, 2)));
    auto mu = (float) ((pow(a, 2) + pow(b, 2)) / (4 * pow(c, 2)));

    vector<std::array<float, 3>> kp = {
            std::array<float, 3>({0.0, 0.0, 0.0}),
            std::array<float, 3>({-mu, mu, static_cast<float>(0.5 - delta)}),
            std::array<float, 3>({mu, -mu, static_cast<float>(0.5 + delta)}),
            std::array<float, 3>({static_cast<float>(0.5 - delta), static_cast<float>(0.5 + delta), -mu}),
            std::array<float, 3>({0.0, 0.5, 0.0}),
            std::array<float, 3>({0.5, 0.0, 0.0}),
            std::array<float, 3>({0.0, 0.0, 0.5}),
            std::array<float, 3>({0.25, 0.25, 0.25}),
            std::array<float, 3>({-zeta, zeta, zeta}),
            std::array<float, 3>({zeta, 1 - zeta, -zeta}),
            std::array<float, 3>({eta, -eta, eta}),
            std::array<float, 3>({1 - eta, eta, -eta}),
            std::array<float, 3>({0.5, 0.5, -0.5}),
    };
    return kp;
}


vector<std::array<float, 3>> orcc(float a, float b, float c) {

    if (c > 0) {}

    auto zeta = float((1 + pow(a, 2) / pow(b, 2)) / 4);
    vector<std::array<float, 3>> kp = {
            std::array<float, 3>({0.0, 0.0, 0.0}),
            std::array<float, 3>({zeta, zeta, 0.5}),
            std::array<float, 3>({-zeta, 1 - zeta, 0.5}),
            std::array<float, 3>({0.0, 0.5, 0.5}),
            std::array<float, 3>({0.0, 0.5, 0.0}),
            std::array<float, 3>({-0.5, 0.5, 0.5}),
            std::array<float, 3>({zeta, zeta, 0.0}),
            std::array<float, 3>({-zeta, 1 - zeta, 0.0}),
            std::array<float, 3>({-0.5, 0.5, 0}),
            std::array<float, 3>({0.0, 0.0, 0.5}),
    };
    return kp;
}


vector<std::array<float, 3>> hex() {


    vector<std::array<float, 3>> kp = {
            std::array<float, 3>({0.0, 0.0, 0.0}),
            std::array<float, 3>({0.0, 0.0, 0.5}),
            std::array<float, 3>({1.0 / 3.0, 1.0 / 3.0, 0.5}),
            std::array<float, 3>({1.0 / 3.0, 1.0 / 3.0, 0.0}),
            std::array<float, 3>({0.5, 0.0, 0.5}),
            std::array<float, 3>({0.5, 0.0, 0.0}),
    };
    return kp;
}


vector<std::array<float, 3>> rhl1(float alpha) {


    float eta = (1 + 4 * cos(alpha)) / (2 + 4 * cos(alpha));
    auto nu = (float) (3.0 / 4.0 - eta / 2.0);

    vector<std::array<float, 3>> kp = {
            std::array<float, 3>({0.0, 0.0, 0.0}),
            std::array<float, 3>({eta, 0.5, static_cast<float>(1 - eta)}),
            std::array<float, 3>({1.0 / 2.0, static_cast<float>(1 - eta), static_cast<float>(eta - 1.0)}),
            std::array<float, 3>({0.5, 0.5, 0.0}),
            std::array<float, 3>({0.5, 0.0, 0.0}),
            std::array<float, 3>({0.0, 0.0, -0.5}),
            std::array<float, 3>({eta, nu, nu}),
            std::array<float, 3>(
                    {static_cast<float>(1.0 - nu), static_cast<float>(1.0 - nu), static_cast<float>(1 - eta)}),
            std::array<float, 3>({nu, nu, static_cast<float>(eta - 1.0)}),
            std::array<float, 3>({static_cast<float>(1.0 - nu), nu, 0.0}),
            std::array<float, 3>({nu, 0.0, -nu}),
            std::array<float, 3>({0.5, 0.5, 0.5}),
    };
    return kp;
}


vector<std::array<float, 3>> rhl2(float alpha) {


    auto eta = float(1 / (2 * pow(tan(alpha / 2.0), 2)));
    auto nu = float(3.0 / 4.0 - eta / 2.0);
    vector<std::array<float, 3>> kp = {
            std::array<float, 3>({0.0, 0.0, 0.0}),
            std::array<float, 3>({0.5, -0.5, 0.0}),
            std::array<float, 3>({0.5, 0.0, 0.0}),
            std::array<float, 3>({1 - nu, -nu, 1 - nu}),
            std::array<float, 3>({nu, static_cast<float>(nu - 1.0), static_cast<float>(nu - 1.0)}),
            std::array<float, 3>({eta, eta, eta}),
            std::array<float, 3>({static_cast<float>(1 - eta), -eta, -eta}),
            std::array<float, 3>({0.5, -0.5, 0.5}),
    };
    return kp;
}


vector<std::array<float, 3>> mcl(float b, float c, float beta) {


    auto eta = float((1 - b * cos(beta) / c) / (2 * pow(sin(beta), 2)));
    auto nu = float(0.5 - eta * c * cos(beta) / b);
    vector<std::array<float, 3>> kp = {
            std::array<float, 3>({0.0, 0.0, 0.0}),
            std::array<float, 3>({0.5, 0.5, 0.0}),
            std::array<float, 3>({0.0, 0.5, 0.5}),
            std::array<float, 3>({0.5, 0.0, 0.5}),
            std::array<float, 3>({0.5, 0.5, -0.5}),
            std::array<float, 3>({0.5, 0.5, 0.5}),
            std::array<float, 3>({0.0, eta, static_cast<float>(1.0 - nu)}),
            std::array<float, 3>({0.0, static_cast<float>(1.0 - eta), nu}),
            std::array<float, 3>({0.0, eta, -nu}),
            std::array<float, 3>({0.5, eta, static_cast<float>(1.0 - nu)}),
            std::array<float, 3>({0.5, 1 - eta, nu}),
            std::array<float, 3>({0.5, 1 - eta, nu}),
            std::array<float, 3>({0.0, 0.5, 0.0}),
            std::array<float, 3>({0.0, 0.0, 0.5}),
            std::array<float, 3>({0.0, 0.0, -0.5}),
            std::array<float, 3>({0.5, 0.0, 0.0}),
    };
    return kp;
}


vector<std::array<float, 3>> mclc1(float a, float b, float c, float alpha) {


    auto zeta = float((2 - b * cos(alpha) / c) / (4 * pow(sin(alpha), 2)));
    auto eta = float(0.5 + 2 * zeta * c * cos(alpha) / b);
    auto psi = float(0.75 - pow(a, 2) / (4 * pow(b, 2) * pow(sin(alpha), 2)));
    auto phi = float(psi + (0.75 - psi) * b * cos(alpha) / c);

    vector<std::array<float, 3>> kp = {
            std::array<float, 3>({0.0, 0.0, 0.0}),
            std::array<float, 3>({0.5, 0.0, 0.0}),
            std::array<float, 3>({0.0, -0.5, 0.0}),
            std::array<float, 3>({1 - zeta, 1 - zeta, 1 - eta}),
            std::array<float, 3>({zeta, zeta, eta}),
            std::array<float, 3>({-zeta, -zeta, 1 - eta}),
            std::array<float, 3>({phi, 1 - phi, 0.5}),
            std::array<float, 3>({1 - phi, phi - 1, 0.5}),
            std::array<float, 3>({0.5, 0.5, 0.5}),
            std::array<float, 3>({0.5, 0.0, 0.5}),
            std::array<float, 3>({1 - psi, psi - 1, 0.0}),
            std::array<float, 3>({psi, 1 - psi, 0.0}),
            std::array<float, 3>({psi - 1, -psi, 0.0}),
            std::array<float, 3>({0.5, 0.5, 0.0}),
            std::array<float, 3>({-0.5, -0.5, 0.0}),
            std::array<float, 3>({0.0, 0.0, 0.5}),
    };
    return kp;
}


vector<std::array<float, 3>> mclc2(float a, float b, float c, float alpha) {


    auto zeta = float((2 - b * cos(alpha) / c) / (4 * pow(sin(alpha), 2)));
    auto eta = float(0.5 + 2 * zeta * c * cos(alpha) / b);
    auto psi = float(0.75 - pow(a, 2) / (4 * pow(b, 2) * pow(sin(alpha), 2)));
    auto phi = float(psi + (0.75 - psi) * b * cos(alpha) / c);

    vector<std::array<float, 3>> kp = {
            std::array<float, 3>({0.0, 0.0, 0.0}),
            std::array<float, 3>({0.5, 0.0, 0.0}),
            std::array<float, 3>({0.0, -0.5, 0.0}),
            std::array<float, 3>({1 - zeta, 1 - zeta, 1 - eta}),
            std::array<float, 3>({zeta, zeta, eta}),
            std::array<float, 3>({-zeta, -zeta, 1 - eta}),
            std::array<float, 3>({1 - zeta, -zeta, 1 - eta}),
            std::array<float, 3>({phi, 1 - phi, 0.5}),
            std::array<float, 3>({1 - phi, phi - 1, 0.5}),
            std::array<float, 3>({0.5, 0.5, 0.5}),
            std::array<float, 3>({0.5, 0.0, 0.5}),
            std::array<float, 3>({1 - psi, psi - 1, 0.0}),
            std::array<float, 3>({psi, 1 - psi, 0.0}),
            std::array<float, 3>({psi - 1, -psi, 0.0}),
            std::array<float, 3>({0.5, 0.5, 0.0}),
            std::array<float, 3>({-0.5, -0.5, 0.0}),
            std::array<float, 3>({0.0, 0.0, 0.5}),
    };
    return kp;
}


vector<std::array<float, 3>> mclc3(float a, float b, float c, float alpha) {


    auto mu = float((1 + pow(b, 2) / pow(a, 2)) / 4.0);
    auto delta = float(b * c * cos(alpha) / (2 * pow(a, 2)));
    auto zeta = float(mu - 0.25 + (1 - b * cos(alpha) / c) / (4 * pow(sin(alpha), 2)));
    auto eta = float(0.5 + 2 * zeta * c * cos(alpha) / b);
    float phi = 1 + zeta - 2 * mu;
    float psi = eta - 2 * delta;
    vector<std::array<float, 3>> kp = {
            std::array<float, 3>({0.0, 0.0, 0.0}),
            std::array<float, 3>({1 - phi, 1 - phi, 1 - psi}),
            std::array<float, 3>({phi, phi - 1, psi}),
            std::array<float, 3>({1 - phi, -phi, 1 - psi}),
            std::array<float, 3>({zeta, zeta, eta}),
            std::array<float, 3>({1 - zeta, -zeta, 1 - eta}),
            std::array<float, 3>({-zeta, -zeta, 1 - eta}),
            std::array<float, 3>({0.5, -0.5, 0.5}),
            std::array<float, 3>({0.5, 0.0, 0.5}),
            std::array<float, 3>({0.5, 0.0, 0.0}),
            std::array<float, 3>({0.0, -0.5, 0.0}),
            std::array<float, 3>({0.5, -0.5, 0.0}),
            std::array<float, 3>({mu, mu, delta}),
            std::array<float, 3>({1 - mu, -mu, -delta}),
            std::array<float, 3>({-mu, -mu, -delta}),
            std::array<float, 3>({mu, mu - 1, delta}),
            std::array<float, 3>({0.0, 0.0, 0.5}),
    };
    return kp;
}


vector<std::array<float, 3>> mclc4(float a, float b, float c, float alpha) {

    auto mu = float((1 + pow(b, 2) / pow(a, 2)) / 4.0);
    auto delta = float(b * c * cos(alpha) / (2 * pow(a, 2)));
    auto zeta = float(mu - 0.25 + (1 - b * cos(alpha) / c) / (4 * pow(sin(alpha), 2)));
    auto eta = float(0.5 + 2 * zeta * c * cos(alpha) / b);
    float phi = 1 + zeta - 2 * mu;
    float psi = eta - 2 * delta;


    vector<std::array<float, 3>> kp = {
            std::array<float, 3>({0.0, 0.0, 0.0}),
            std::array<float, 3>({1 - phi, 1 - phi, 1 - psi}),
            std::array<float, 3>({phi, phi - 1, psi}),
            std::array<float, 3>({1 - phi, -phi, 1 - psi}),
            std::array<float, 3>({zeta, zeta, eta}),
            std::array<float, 3>({1 - zeta, -zeta, 1 - eta}),
            std::array<float, 3>({-zeta, -zeta, 1 - eta}),
            std::array<float, 3>({0.5, -0.5, 0.5}),
            std::array<float, 3>({0.5, 0.0, 0.5}),
            std::array<float, 3>({0.5, 0.0, 0.0}),
            std::array<float, 3>({0.0, -0.5, 0.0}),
            std::array<float, 3>({0.5, -0.5, 0.0}),
            std::array<float, 3>({mu, mu, delta}),
            std::array<float, 3>({1 - mu, -mu, -delta}),
            std::array<float, 3>({-mu, -mu, -delta}),
            std::array<float, 3>({mu, mu - 1, delta}),
            std::array<float, 3>({0.0, 0.0, 0.5}),
    };
    return kp;
}


vector<std::array<float, 3>> mclc5(float a, float b, float c, float alpha) {


    auto zeta = float((pow(b, 2) / pow(a, 2) + (1 - b * cos(alpha) / c) / pow(sin(alpha), 2)) / 4);
    auto eta = float(0.5 + 2 * zeta * c * cos(alpha) / b);
    auto mu = float(eta / 2 + pow(b, 2) / (4 * pow(a, 2)) - b * c * cos(alpha) / (2 * pow(a, 2)));
    float nu = 2 * mu - zeta;
    auto rho = float(1 - zeta * pow(a, 2) / pow(b, 2));
    auto omega = float((4 * nu - 1 - pow(b, 2) * pow(sin(alpha), 2) / pow(a, 2)) * c / (2 * b * cos(alpha)));
    auto delta = float(zeta * c * cos(alpha) / b + omega / 2 - 0.25);

    vector<std::array<float, 3>> kp = {
            std::array<float, 3>({0.0, 0.0, 0.0}),
            std::array<float, 3>({nu, nu, omega}),
            std::array<float, 3>({1 - nu, 1 - nu, 1 - omega}),
            std::array<float, 3>({nu, nu - 1, omega}),
            std::array<float, 3>({zeta, zeta, eta}),
            std::array<float, 3>({1 - zeta, -zeta, 1 - eta}),
            std::array<float, 3>({-zeta, -zeta, 1 - eta}),
            std::array<float, 3>({rho, 1 - rho, 0.5}),
            std::array<float, 3>({1 - rho, rho - 1, 0.5}),
            std::array<float, 3>({0.5, 0.5, 0.5}),
            std::array<float, 3>({0.5, 0.0, 0.5}),
            std::array<float, 3>({0.5, 0.0, 0.0}),
            std::array<float, 3>({0.0, -0.5, 0.0}),
            std::array<float, 3>({0.5, -0.5, 0.0}),
            std::array<float, 3>({mu, mu, delta}),
            std::array<float, 3>({1 - mu, -mu, -delta}),
            std::array<float, 3>({-mu, -mu, -delta}),
            std::array<float, 3>({mu, mu - 1, delta}),
            std::array<float, 3>({0.0, 0.0, 0.5}),
    };
    return kp;
}


vector<std::array<float, 3>> tria() {


    vector<std::array<float, 3>> kp = {
            std::array<float, 3>({0.0, 0.0, 0.0}),
            std::array<float, 3>({0.5, 0.5, 0.0}),
            std::array<float, 3>({0.0, 0.5, 0.5}),
            std::array<float, 3>({0.5, 0.0, 0.5}),
            std::array<float, 3>({0.5, 0.5, 0.5}),
            std::array<float, 3>({0.5, 0.0, 0.0}),
            std::array<float, 3>({0.0, 0.5, 0.0}),
            std::array<float, 3>({0.0, 0.0, 0.5}),
    };
    return kp;
}


vector<std::array<float, 3>> trib() {


    vector<std::array<float, 3>> kp = {
            std::array<float, 3>({0.0, 0.0, 0.0}),
            std::array<float, 3>({0.5, -0.5, 0.0}),
            std::array<float, 3>({0.0, 0.0, 0.5}),
            std::array<float, 3>({-0.5, -0.5, 0.5}),
            std::array<float, 3>({0.0, -0.5, 0.5}),
            std::array<float, 3>({0.0, -0.5, 0.0}),
            std::array<float, 3>({0.5, 0.0, 0.0}),
            std::array<float, 3>({-0.5, 0.0, 0.5}),
    };
    return kp;
}

//default zero
vector<std::array<float, 3>> default_zeros() {


    vector<std::array<float, 3>> kp = {
            std::array<float, 3>({0.0, 0.0, 0.0}),
            std::array<float, 3>({0.0, -0.0, 0.0}),

    };
    return kp;
}


std::vector<std::array<float, 3>> kp(float ty, float a, float b, float c, float alpha_prim, float alpha_conv) {
    int t = (int) ty;
    float pi = 3.1415926;
    switch (t) {
        case 0:
            return cubic();
        case 1:
            return fcc();
        case 2:
            return bcc();
        case 3:
            return tet();
        case 4 :
            return bctet1(c, a); // 4  conv.lattice.abc
        case 5 :
            return bctet2(c, a); // 5  conv.lattice.abc
        case 6 :
            return orc(); // 6 /
        case 7:
            return orcf1(a, b, c); // 7
        case 8 :
            return orcf2(a, b, c); // 8
        case 9:
            return orcf3(a, b, c); // 9
        case 10 :
            return orci(a, b, c); // 10
        case 11:
            return orcc(a, b, c); // 11
        case 12 :
            return hex(); // 12 /
        case 13 :
            return rhl1(alpha_prim * pi / 180); // 13  _prim.lattice.parameters
        case 14:
            return rhl2(alpha_prim * pi / 180); // 14 _prim.lattice.parameters
        case 15:
            return mcl(b, c, alpha_conv * pi / 180); // 15  _conv.lattice.parameters
        case 16:
            return mclc1(a, b, c, alpha_conv * pi / 180); // 16
        case 17:
            return mclc2(a, b, c, alpha_conv * pi / 180); // 17
        case 18:
            return mclc3(a, b, c, alpha_conv * pi / 180); // 18
        case 19:
            return mclc4(a, b, c, alpha_conv * pi / 180); // 19
        case 20:
            return mclc5(a, b, c, alpha_conv * pi / 180); // 20
        case 21:
            return tria(); // 21 /
        case 22:
            return trib(); // 22 /
        case 23:
            return tria(); // 23 /
        case 24:
            return trib(); // 24 /

        default:
            return default_zeros();

    }
}

std::array<std::array<float, 3>, 20> kpuni(float ty, float a, float b, float c, float alpha_prim, float alpha_conv) {

    std::array<std::array<float, 3>, 20> temp = {std::array<float, 3>({0.0, 0.0, 0.0}),};
    auto kps = kp(ty, a, b, c, alpha_prim, alpha_conv);
    std::copy(std::begin(kps), std::end(kps), std::begin(temp));
    return temp;
}

typedef Eigen::Matrix<float, 3, 20> Mat20;
typedef std::vector<Mat20> Mats20;
typedef Eigen::Matrix<float, 1, 6> Vec6;
typedef Eigen::Matrix<float, Eigen::Dynamic, 6> Vecs6;

Mat20 kpeuni(Vec6 arr) {

    Mat20 temp = Mat20::Zero();
    auto kps = kp(arr(0, 0), arr(0, 1), arr(0, 2), arr(0, 3), arr(0, 4), arr(0, 5));

    for (int i = 0; i < kps.size(); i++) {
        temp.col(i) = Map<Eigen::Matrix<float, 1, 3>>(&kps[i][0], kps[i].size());
    }
    return temp;
}


Mats20 kps(Vecs6 arr) {
    Mats20 ks;

    for (int i = 0; i < arr.rows(); i++) {
        auto si = arr.row(i);
        Mat20 k = kpeuni(si);
        ks.push_back(k);
    }
    return ks;

}


PYBIND11_MODULE(kpcpu, m) {
    m.doc() = "A function which return high symmetry K Points of cpu"; // optional module docstring

    m.def("kp", &kp, "A function which return high symmetry K Points.");
    m.def("kpuni", &kpuni, "A function which return high symmetry K Points with the size 20x3,(list type).");
    m.def("kpeuni", &kpeuni, "A function which return high symmetry K Points with the size 20x3(np.ndarray type).");
    m.def("kps", &kps, "A function which return high symmetry K Points with the size c*[20x3] (list of np.ndarray type).");
}
