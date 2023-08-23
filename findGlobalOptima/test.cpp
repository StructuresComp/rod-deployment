#include <iostream>
#include <fstream>
#include <string>
#include "eigenIncludes.h"


int main()
{
	ifstream myfile("simDER.txt"); 

	int row1 =  101; 

	MatrixXd data = MatrixXd(row1, 4);
    double a ;
	if (myfile.is_open())
	{
		for (int i = 0; i<row1* 4; i++)
		{
			myfile>>a;
			// cout<<a<<endl;
			if (i%4 == 0)
				data(i/4, 0) = a;
			else if (i%4 == 1)
				data(i/4, 1) = a;
			else if (i%4 == 2)
				data(i/4, 2) = a;
			else if (i%4 == 3)
				data(i/4, 3) = a;
		}
	}

	cout<< data.block<<endl;

	return 0;

}