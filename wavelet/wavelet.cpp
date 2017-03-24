#include "wavelet.h"
using namespace std;
void wavelet::wavedec(double* srcdata, int dataLen, double* det, double* app)
{
	int filterLen = 4; // �˲�������


					   //  д���˲���
	double filter_LD[4] = { -0.12941, 0.22414, 0.83652, 0.48296 };// �ֽ� ��ͨ���� DET
	double filter_HD[4] = { -0.48296,0.83652,-0.22414,-0.12941 }; // �ֽ� ��ͨ���� APP


	int exLen = dataLen + 2 * filterLen - 2;  // ���غ�ĳ���
											  // �Գ����غ������
	double* exdata = new double[exLen];   //  ��ʱ���غ������
										  //  ���Գ�����
	for (int i = 0; i < exLen; i++)
	{
		if (i < filterLen - 1)
			exdata[i] = srcdata[filterLen - 1 - i - 1];
		else if (i >= filterLen - 1 && i < exLen - filterLen + 1)
			exdata[i] = srcdata[i - 3];
		else
			exdata[i] = srcdata[dataLen - (i - dataLen - filterLen + 1) - 1];

	}


	//  �����
	int convLen = exLen - filterLen + 1; // ����ÿ��������������ĳ���
	double* pre_det = new double[convLen];  // ������ϸ�ڣ���ͨ������
	double* pre_app = new double[convLen];  // ��������������ͨ������
	for (int i = 0; i < convLen; i++)
	{
		// ������ݳ�ʼ��
		pre_det[i] = 0.0;  //  ��ͨ����
		pre_app[i] = 0.0;  //  ��ͨ����
		for (int j = filterLen - 1, k = 0; j >= 0; j--, k++)
		{
			pre_det[i] += filter_HD[j] * exdata[i + k];  //   ���ͨ�˲��������
			pre_app[i] += filter_LD[j] * exdata[i + k];  //   ���ͨ�˲��������
		}
	}
	delete[] exdata;    //  �ͷ���ʱ���غ������


						//  ���²���
	int detLen = convLen / 2;   //  �������ϸ�ڳ���
	for (int i = 1, j = 0; i < convLen; i += 2, j++)
	{
		det[j] = pre_det[i];
		app[j] = pre_app[i];
	}
	//------------------�������----------------
	/*
	cout << detLen << endl;
	for (int i = 0; i < detLen; i++)
	{
	cout << det[i] << "  " << app[i] << endl;
	}
	*/
	delete[] pre_app;   //  �ͷſռ�
	delete[] pre_det;

}

void wavelet::waverec(double* det, double* app, int n, double* outdata)
{
	int filterLen = 4; // �˲�������

					   //  д���˲���
	double filter_HR[4] = { -0.12941,-0.22414,0.83652,-0.48296 }; // �ع� ��ͨ���� DET
	double filter_LR[4] = { 0.48296,0.83652,0.22414,-0.12941 };   // �ع� ��ͨ���� APP


																  //  ׼���ع�
	int upLen = n * 2 + 1;  //  �ϲ�����ĳ���
	double* up_det = new double[upLen];
	double* up_app = new double[upLen];
	//  �ϲ���
	for (int i = 0; i < n; i++)
	{
		up_det[2 * i] = 0.0;
		up_det[2 * i + 1] = det[i];
		up_app[2 * i] = 0.0;
		up_app[2 * i + 1] = app[i];
	}
	up_det[2 * n] = up_app[2 * n] = 0.0; //  ĩλ����
										 //-----------------�������--------------
										 /*cout << "up_det: ";
										 for (int i = 0; i < upLen; i++)
										 {
										 cout << up_det[i]<<" ";
										 }
										 cout << endl;

										 */
										 //  ����
	int exLen = upLen + 2 * filterLen - 2;  //  ���غ�����鳤��
	double* exdet = new double[exLen];   //  ���غ��det����
	double* exapp = new double[exLen];   //  ���غ��app����
	for (int i = 0; i < exLen; i++)
	{
		if (i < filterLen - 1)
		{
			exdet[i] = up_det[filterLen - 1 - i - 1];
			exapp[i] = up_app[filterLen - 1 - i - 1];
		}

		else if (i >= filterLen - 1 && i < exLen - filterLen + 1)
		{
			exdet[i] = up_det[i - 3];
			exapp[i] = up_app[i - 3];
		}

		else
		{
			exdet[i] = up_det[upLen - (i - upLen - filterLen + 1) - 1];
			exapp[i] = up_app[upLen - (i - upLen - filterLen + 1) - 1];
		}
	}
	//-----------------�������--------------
	/*cout << "exdet: ";
	for (int i = 0; i < exLen; i++)
	{
	cout << exdet[i] << " ";
	}
	cout << endl;
	*/
	delete[] up_det;
	delete[] up_app;
	//  �����
	int aftLen = exLen - filterLen + 1;  //   �����ĳ���
	double* aft_det = new double[aftLen];     //  ������det����
	double* aft_app = new double[aftLen];     //  ������app����
											  //double* recdata = new double[aftLen];   //   �ع��������
	for (int i = filterLen - 1, j = 0; i < aftLen - filterLen + 1; i++)
	{
		// ������ݳ�ʼ��
		aft_det[i] = 0.0;  //  ��ͨ����
		aft_app[i] = 0.0;  //  ��ͨ����
		for (int j = filterLen - 1, k = 0; j >= 0; j--, k++)
		{
			aft_det[i] += filter_HR[j] * exdet[i + k];  //   ���ͨ�˲��������
			aft_app[i] += filter_LR[j] * exapp[i + k];  //   ���ͨ�˲��������
		}
		//cout << aft_det[i] << " " << aft_app[i] << endl;
		outdata[j++] = aft_det[i] + aft_app[i];  //  ��ͨ���ͨ��������ع�
	}
	delete[] exdet;
	delete[] exapp;
	delete[] aft_det;
	delete[] aft_app;

}

double wavelet::getT(double* det, int detLen)
{
	//  ��ȡ��ֵ
	double sigma = 0.0;  //  ��ʼ��sigmaֵ
						 //-----------------�������--------------
						 //cout << "abs(det): ";
						 //for (int i = 0; i < detLen; i++)
						 //{
						 //	if (det[i] < 0)
						 //		det[i] = -det[i];
						 //	cout << det[i]<<" ";
						 //}
						 //cout << endl;
	double* temp = new double[detLen];
	for (int i = 0; i < detLen; i++)
	{
		if (det[i] < 0)
			temp[i] = -det[i];
		else
			temp[i] = det[i];
	}
	sort(temp, temp + detLen);
	sigma = temp[detLen / 2];
	cout << "sigma: " << sigma << endl;
	//  ������ֵ T
	double T = sigma*sqrt(2 * log(detLen));
	cout << "T: " << T << endl;
	return T;
}

void wavelet::wavecle(double* det, int detLen, double T)
{
	//double T = getT(det, detLen);
	//  ʹ����ֵ����
	for (int i = 0; i < detLen; i++)
	{
		if (det[i] < T)
			det[i] = 0.0;
	}
}

void wavelet::wavecle(double* srcdata, int dataLen, int num, double* outdata)
{
	vector<double*> app;
	vector<double*> det;
	int length = dataLen;
	double* data = srcdata;

	//  �ֽⲿ��
	for (int i = 0; i<num; i++)
	{
		det[i] = new double[(length + 4 - 1) / 2];
		app[i] = new double[(length + 4 - 1) / 2];
		wavedec(data, dataLen, det[i], app[i]);
		data = app[i];
		length = (length + 4 - 1) / 2;
	}
	//  ��ֵȥ��
	int T = getT(det[0], (dataLen + 4 - 1) / 2);
	int detlen = (dataLen + 4 - 1) / 2;
	for (int i = 0; i < num; i++)
	{
		wavecle(det[i], detlen, T);
		detlen = (detlen + 4 - 1) / 2;
	}
	//  �ع�����
	for (int i = num - 1; i >= 0; i--)
	{
		waverec(det[i], app[i], detlen, outdata);
		if (i != 0)
			app[i] = outdata;
		detlen = 2 * detlen - 4 + 1;
	}

}

void exchange(double &a, double &b)
{
	double temp;
	temp = a;
	a = b;
	b = temp;
}