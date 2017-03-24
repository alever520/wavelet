#pragma once
#pragma once
#include<iostream>
#include<algorithm>
#include<vector>
void exchange(double &a, double &b);
class wavelet
{

public:
	//wavelet(double* in, int n);
	void wavedec(double* srcdata, int dataLen, double* det, double* app);
	void waverec(double* det, double* app, int n, double* outdata);
	void wavecle(double* det, int n, double T);
	double getT(double* det, int detLen);
	void wavecle(double* srcdata, int dataLen, int num, double* outdata);

};