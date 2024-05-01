
#include"mpi.h"
#include <vector>
#include<iostream>
#include<fstream>
#include<cmath>
#include<omp.h>
#include<iomanip>
#include<cmath>

const double PI = 3.141592653589793;
double Get_New_Point_In_Time(double prev_point, double curr_point, double next_point, double deltaT, double deltaX, double sigma)
{
  double tmp = 0.0;
  tmp = prev_point - 2 * curr_point + next_point;

  return((sigma * deltaT) / (deltaX * deltaX) * tmp + curr_point);
}

void Insert_In_File(double* new_array, size_t count, std::fstream& stream)
{
  if (stream.is_open())
  {
    for (size_t i = 0; i < count; i++)
      stream << new_array[i] << std::setprecision(15) << " ";
    stream << "\n";
    stream << "\n";
  }
  else
    throw "The file isn't exist";
}

int main(int argc, char** argv) {

 
  int myrank;
  int count_of_proc;
  long int count = 0;
  long int points_per_proc = 0;
  int remainder = 0;
  int q = 0;
  double time = 0.0;
  double X_max = 0.0;
  double deltaX = 0.0;
  double T_max = 0.0;
  double deltaT = 0.0;
  double sigma = 0.0;
  double epsilon = 0.000001;
  size_t iteration = 0;

  std::fstream f;
  f.open("output.txt", std::fstream::in | std::fstream::out);
  std::ifstream f_in;
  //////////////////////////////////////////////////////////////////////////////////////////
  f_in.open("const_initial.txt");
  if (f_in.is_open())
  {
    f_in >> X_max;
    f_in >> deltaX;
    f_in >> T_max;
    f_in >> deltaT;
    f_in >> sigma;
    f_in >> epsilon;

  }
  else
    throw "file wasn't opened";
  f_in.close();


  MPI_Init(&argc, &argv);
  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &count_of_proc);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  double* array = nullptr;
  double* get_array = nullptr;
  double* temp_result_array = nullptr;

  int size_byte = 0;
  double left = 0, right = 0;

  std::vector<int> A1(count_of_proc); // size_of scatterv
  std::vector<int> A2(count_of_proc); // size_of gatherv
  std::vector<int> B1(count_of_proc); // displs for scatterv
  std::vector<int> B2(count_of_proc); //displs for gatherv

  count = (X_max / deltaX) + 1;

  points_per_proc = count / count_of_proc;
  remainder = count % count_of_proc;

  if (myrank == count_of_proc - 1)
    points_per_proc = points_per_proc + remainder;

  temp_result_array = new double[points_per_proc];
  std::cout <<"rank = "<<myrank << "points_per_proc" << points_per_proc << std::endl;
  if (myrank == 0)
  {
    array = new double[count];
    for (int i = 0; i < count; i++)
    {
     array[i] = sin(PI * deltaX * i);
     // array[i] = i;
      std::cout << array[i] << ' ';
    }
    array[count - 1] = 0.0;
    std::cout << std::endl;
    Insert_In_File(array, count, f);
    for (int q = 0; q < count_of_proc; q++)
    {
      A1[q] = points_per_proc;
      A2[q] = points_per_proc;
      B1[q] = 0;
      B2[q] = 0;
    }
    A1[count_of_proc - 1] += remainder;
    A2[count_of_proc - 1] += remainder;
    for (int q = 1; q < count_of_proc; q++)
    {
      B1[q] = B1[q - 1] + A1[q - 1];
      B2[q] = B2[q - 1] + A2[q - 1];
    }

    get_array = new double[A2[0]];
  }
  else
  {
      get_array = new double[points_per_proc];
  }


  while (time < T_max)
  {
    
    MPI_Scatterv(array, A1.data(), B1.data(), MPI_DOUBLE, get_array, points_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (myrank == 0)
    {
      MPI_Send(array + points_per_proc - 1, 1, MPI_DOUBLE, myrank + 1, myrank, MPI_COMM_WORLD);
      MPI_Recv(&right, 1, MPI_DOUBLE, myrank + 1, myrank + 1, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
    }
    else if (myrank == (count_of_proc - 1))
    {
      MPI_Recv(&left, 1, MPI_DOUBLE, myrank - 1, myrank - 1, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
      MPI_Send(get_array, 1, MPI_DOUBLE, myrank - 1, myrank, MPI_COMM_WORLD);

    }
    else
    {
      MPI_Recv(&left, 1, MPI_DOUBLE, myrank - 1, myrank - 1, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
      MPI_Send(get_array + points_per_proc - 1, 1, MPI_DOUBLE, myrank + 1, myrank, MPI_COMM_WORLD);
      MPI_Recv(&right, 1, MPI_DOUBLE, myrank + 1, myrank + 1, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
      MPI_Send(get_array, 1, MPI_DOUBLE, myrank - 1, myrank, MPI_COMM_WORLD);
    }
   // std::cout << myrank << "|" << " " << left << ' ' << right <<" iteration = "<<iteration << std::endl;

     /* for (int q = 0; q < points_per_proc; q++)
        std::cout << "my rank = " << myrank << ' ' << get_array[q] << " iteration = " << iteration << ' ';
      std::cout<<std::endl;*/
    
      iteration++;
    if (myrank == 0)
    {
      
      temp_result_array[0] = 0.0;
      for (size_t q = 1; q < points_per_proc - 1; q++)
      {
        temp_result_array[q] = Get_New_Point_In_Time(get_array[q - 1], get_array[q], get_array[q + 1], deltaT, deltaX, sigma);
        // std::cout << get_array[q - 1] << ' ' << get_array[q] << ' ' << get_array[q + 1] << std::endl;
      }
      // std::cout << get_array[points_per_proc - 2] << ' ' << get_array[points_per_proc-1] << ' ' << right << std::endl;
      temp_result_array[points_per_proc - 1] = Get_New_Point_In_Time(get_array[points_per_proc - 2], get_array[points_per_proc - 1], right, deltaT, deltaX, sigma);

    }

    else if (myrank == count_of_proc - 1)
    {
      temp_result_array[points_per_proc - 1] = 0;
      temp_result_array[0] = Get_New_Point_In_Time(left, get_array[0], get_array[1], deltaT, deltaX, sigma);
      for (size_t q = 1; q < points_per_proc - 1; q++)
      {
        temp_result_array[q] = Get_New_Point_In_Time(get_array[q - 1], get_array[q], get_array[q + 1], deltaT, deltaX, sigma);
      }
     
    }

    else
    {
    
      temp_result_array[0] = Get_New_Point_In_Time(left, get_array[0], get_array[ 1], deltaT, deltaX, sigma);
      temp_result_array[points_per_proc-1] = Get_New_Point_In_Time(get_array[points_per_proc-2], get_array[points_per_proc-1], right, deltaT, deltaX, sigma);
      
      for (size_t q = 1; q < points_per_proc - 1; q++)
      {
        temp_result_array[q] = Get_New_Point_In_Time(get_array[q - 1], get_array[q], get_array[q + 1], deltaT, deltaX, sigma);
       
      }
     
    }
  
    MPI_Gatherv(temp_result_array, points_per_proc, MPI_DOUBLE, array, A2.data(), B2.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    

    time += deltaT;
 
  } 
 

  //f.close();
  if (myrank == 0)
  {
    for (int q = 0; q < count; q++)
      std::cout << array[q] << ' ';
    std::cout << std::endl;
    Insert_In_File(array, count, f);
  }
  std::cout << "myrank = " << myrank << std::endl;
  delete[] array;
  delete[] temp_result_array;
  delete[] get_array;
  MPI_Finalize();
 

  return 0;
}