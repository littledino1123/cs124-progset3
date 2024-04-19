#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <queue>
#include <chrono>
#include <random>
#include <climits>

using namespace std;
using namespace std::chrono;
mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
uniform_real_distribution<float> distribution(0.0, 1.0);

// set Max iter to 25000
static const int MAX_ITER = 25000;

// Define Karmarkar-Karp function
long long karmarkarKarp(const std::vector<long long> &numbers);

// Define repeated random function
long long repeatedRandom(const std::vector<long long> &numbers, bool prepartition);

// Define hill climbing function
long long hillClimbing(const std::vector<long long> &numbers, bool prepartition);

// Define sim annealing function
long long simulatedAnnealing(const std::vector<long long> &numbers, bool prepartition);

// T function for simulated annealing
long long T(int iter) {
  return pow(10,10) * pow(0.8, iter / 300);
}

int main(int argc, char *argv[])
{
  if (argc != 4)
  {
    std::cerr << "Usage: " << argv[0] << " flag algorithm inputfile" << std::endl;
    return 1;
  }

  int flag = atoi(argv[1]);
  // case for generating 50 random instances
  if (flag == 4) {
    long double kkTime = 0;
    long double rrTime = 0;
    long double hcTime = 0;
    long double saTime = 0;
    long double rrPTime = 0;
    long double hcPTime = 0;
    long double saPTime = 0;

    for (int i = 0; i < 50; ++i) {
      vector<long long> newNumbers;
      for (int j = 0; j < 100; ++j) {
        long long random_number = ((long long)rand() << 30) | rand(); // Shift by 30 bits to form a 60-bit integer
        long long scaled_random_number = random_number % 1000000000000LL; // Ensure the number is within [0, 10^12 - 1]
        newNumbers.push_back(scaled_random_number + 1);
      }
      auto start = high_resolution_clock::now();
      cout << "KK (no partition)" << karmarkarKarp(newNumbers) << endl;
      auto stop = high_resolution_clock::now();
      auto duration = duration_cast<milliseconds>(stop - start);
      kkTime += duration.count() / 1000.0;

      start = high_resolution_clock::now();
      cout << "RR (no partition)" << repeatedRandom(newNumbers, false) << endl;
      stop = high_resolution_clock::now();
      duration = duration_cast<milliseconds>(stop - start);
      rrTime += duration.count() / 1000.0;

      start = high_resolution_clock::now();
      cout << "HC (no partition)" << hillClimbing(newNumbers, false) << endl;
      stop = high_resolution_clock::now();
      duration = duration_cast<milliseconds>(stop - start);
      hcTime += duration.count() / 1000.0;

      start = high_resolution_clock::now();
      cout << "SA (no partition)" << simulatedAnnealing(newNumbers, false) << endl;
      stop = high_resolution_clock::now();
      duration = duration_cast<milliseconds>(stop - start);
      saTime += duration.count() / 1000.0;

      start = high_resolution_clock::now();
      cout << "RR (partition)" << repeatedRandom(newNumbers, true) << endl;
      stop = high_resolution_clock::now();
      duration = duration_cast<milliseconds>(stop - start);
      rrPTime += duration.count() / 1000.0;
      
      start = high_resolution_clock::now();
      cout << "HC (partition)" << hillClimbing(newNumbers, true) << endl;
      stop = high_resolution_clock::now();
      duration = duration_cast<milliseconds>(stop - start);
      hcPTime += duration.count() / 1000.0;

      start = high_resolution_clock::now();
      cout << "SA (partition)" << simulatedAnnealing(newNumbers, true) << endl;
      stop = high_resolution_clock::now();
      duration = duration_cast<milliseconds>(stop - start);
      saPTime += duration.count() / 1000.0;
    }
    cout << "KK (no partition) average time: " << kkTime / 50 << endl;
    cout << "RR (no partition) average time: " << rrTime / 50 << endl;
    cout << "HC (no partition) average time: " << hcTime / 50 << endl;
    cout << "SA (no partition) average time: " << saTime / 50 << endl;
    cout << "RR (partition) average time: " << rrPTime / 50 << endl;
    cout << "HC (partition) average time: " << hcPTime / 50 << endl;
    cout << "SA (partition) average time: " << saPTime / 50 << endl;
    return 0;
  }
  int algorithm = atoi(argv[2]);
  string inputfile = argv[3];

  // Read the integers from the input file.
  ifstream infile(inputfile);
  if (!infile)
  {
    cerr << "Error opening file: " << inputfile << endl;
    return 1;
  }

  std::vector<long long> numbers;
  long long number;
  while (infile >> number)
  {
    numbers.push_back(number);
  }

  // breaks down into cases
  long long residue;
  switch (algorithm)
  {
  case 0:
    residue = karmarkarKarp(numbers);
    break;
  case 1:
    residue = repeatedRandom(numbers, false);
    break;
  case 2:
    residue = hillClimbing(numbers, false);
    break;
  case 3:
    residue = simulatedAnnealing(numbers, false);
    break;
  case 11:
    residue = repeatedRandom(numbers, true);
    break;
  case 12:
    residue = hillClimbing(numbers, true);
    break;
  case 13:
    residue = simulatedAnnealing(numbers, true);
    break;
  
  default:
    std::cerr << "Invalid algorithm code." << std::endl;
    return 1;
  }

  std::cout << residue << std::endl;
  return 0;
}

long long karmarkarKarp(const std::vector<long long> &numbers)
{
  std::priority_queue<long long> heap;
  for (auto number : numbers)
  {
    heap.push(number);
  }

  while (heap.size() > 1)
  {
    long long largest = heap.top();
    heap.pop();
    long long second_largest = heap.top();
    heap.pop();
    heap.push(largest - second_largest);
  }

  return heap.empty() ? 0 : heap.top();
}

long long repeatedRandom(const std::vector<long long> &numbers, bool prepartition) {
  long long bestResidue = INT_MAX;
  if (prepartition){
      // Generate a random solution 
      std::vector<long long> partition(numbers.size());
      for (int i = 0; i < numbers.size(); ++i) {
        partition[i] = rand() % numbers.size() + 1;
      }
      // Calculate the first residue 
      std::vector<long long> newNumbers(numbers.size(), 0);
      for (int i = 0; i < numbers.size(); ++i) {
        newNumbers[partition[i] - 1] += numbers[i];
      }
      long long newResidue = karmarkarKarp(newNumbers);
      bestResidue = std::min(abs(bestResidue), abs(newResidue));
      
      // loop through MAX_ITER times
      for (int iter = 0; iter < MAX_ITER; ++iter){
        std::vector<long long> newPartition(numbers.size());
        for (int i = 0; i < numbers.size(); ++i) {
          newPartition[i] = rand() % numbers.size() + 1;
        }
        std::vector<long long> newNumbers(numbers.size(), 0);
        for (int i = 0; i < numbers.size(); ++i) {
          newNumbers[newPartition[i] - 1] += numbers[i];
        }
        long long newResidue = karmarkarKarp(newNumbers);
        bestResidue = std::min(abs(bestResidue), abs(newResidue));
      }
  }
  // No prepartition
  else{
     
    for (int iter = 0; iter < MAX_ITER; ++iter) {
      std::vector<long long> newNumbers = numbers;

      // Generate a random solution
      for (auto &number : newNumbers){
        float w = distribution(rng);
        if (w > 0.5) {
          number = -number;
        }
      }
      long long newResidue = 0;
      for (auto number : newNumbers) {
              newResidue += number; 
          }
      bestResidue = std::min(abs(bestResidue), abs(newResidue));
    }

  
  }
  return bestResidue;
  
}

long long hillClimbing(const std::vector<long long> &numbers, bool prepartition)
{
  long long bestResidue = INT_MAX;
  if (prepartition){
    // Generate a random solution 
    std::vector<long long> partition(numbers.size());
    for (int i = 0; i < numbers.size(); ++i) {
      partition[i] = rand() % numbers.size() + 1;
    }

    // Calculate the first residue 
    std::vector<long long> newNumbers(numbers.size(), 0);
    for (int i = 0; i < numbers.size(); ++i) {
      newNumbers[partition[i] - 1] += numbers[i]; 
    }
    bestResidue = karmarkarKarp(newNumbers);
    
    for (int iter = 0; iter < MAX_ITER; ++iter){
      std::vector<long long> newPartition = partition;
      // find neighbors
      int i = rand() % newPartition.size();
      int j = rand() % newPartition.size() + 1;
      // make sure that partition[i] != j
      while (partition[i] == j) 
      {
        j = rand() % newPartition.size() + 1;
      }
      newPartition[i] = j;
      
      std::vector<long long> newNumbers(numbers.size(), 0);
      for (int i = 0; i < numbers.size(); ++i) {
        newNumbers[newPartition[i] - 1] += numbers[i]; 
      }
      long long newResidue = karmarkarKarp(newNumbers);
      bestResidue = std::min(abs(bestResidue), abs(newResidue));
    }
  }
  else{
    srand(time(NULL));
    std::vector<long long> currentNumbers = numbers;
    for (auto &number : currentNumbers) {
        float w = distribution(rng);
        if (w > 0.5)
        {
          number = -number; 
        }
    }
    bestResidue = 0;
    for (auto number : currentNumbers) {
              bestResidue += number;
    }
    bestResidue = abs(bestResidue);

    for (int iter = 0; iter < MAX_ITER; ++iter) {
      std::vector<long long> newNumbers = currentNumbers;
      int i = rand() % newNumbers.size();
      int j = rand() % newNumbers.size();
      // Ensure that i and j are different 
      while (i == j) {
        j = rand() % newNumbers.size();
      }
      newNumbers[i] *= -1; 
      float w = distribution(rng);
      // flip j with probability 0.5
      if (w > 0.5) {
        newNumbers[j] *= -1; 
      }

      long long newResidue = 0;
      for (auto number : newNumbers) {
              newResidue += number; 
      }
      newResidue = abs(newResidue);
      if (newResidue < bestResidue) {
        bestResidue = newResidue;
        currentNumbers = newNumbers;
      }
    }
  }
  return bestResidue;
}

long long simulatedAnnealing(const std::vector<long long> &numbers, bool prepartition)
{
  long long bestResidue = INT_MAX;
  srand(time(NULL));

  if (prepartition) {
    // Generate a random solution 
    std::vector<long long> partition(numbers.size());
    for (int i = 0; i < numbers.size(); ++i) {
      partition[i] = rand() % numbers.size() + 1;
    }

    // Calculate the first residue 
    std::vector<long long> newNumbers(numbers.size(), 0);
    for (int i = 0; i < numbers.size(); ++i) {
      newNumbers[partition[i] - 1] += numbers[i]; 
    }
    long long currentResidue = karmarkarKarp(newNumbers);
    bestResidue = currentResidue;
    
    for (int iter = 0; iter < MAX_ITER; ++iter) {
      std::vector<long long> newPartition = partition;
      int i = rand() % newPartition.size();
      int j = rand() % newPartition.size() + 1;
      // Ensure that partition[i] != j
      while (partition[i] == j) {
        j = rand() % newPartition.size() + 1;
      }
      newPartition[i] = j; 
      
      std::vector<long long> newNumbers(numbers.size(), 0);
      for (int i = 0; i < numbers.size(); ++i) {
        newNumbers[newPartition[i] - 1] += numbers[i];
      }
      long long newResidue = karmarkarKarp(newNumbers);
      if (newResidue < currentResidue){
        currentResidue = newResidue;
        partition = newPartition;
      }
      else if (exp(-(newResidue - currentResidue) / T(iter)) > distribution(rng)){
        currentResidue = newResidue;
        partition = newPartition;
      }
      if (currentResidue < bestResidue) {
        bestResidue = currentResidue;
        partition = newPartition;
      }
    }

  }
  else {
    
    std::vector<long long> currentNumbers = numbers;
    for (auto &number : currentNumbers) {
        float w = distribution(rng);
        if (w > 0.5) {
          number = -number; 
        }
    }
    long long currentResidue = 0;
    for (auto number : currentNumbers) {
              currentResidue += number; // Add each number to the sum
          }
    currentResidue = abs(currentResidue);

    bestResidue = currentResidue;

    for (int iter = 0; iter < MAX_ITER; ++iter) {
      std::vector<long long> newNumbers = currentNumbers;
      int i = rand() % newNumbers.size();
      int j = rand() % newNumbers.size();
      // Ensure that i and j are different
      while (i == j) 
      {
        j = rand() % newNumbers.size();
      }
      newNumbers[i] *= -1; 
      float w = distribution(rng);
      if (w > 0.5) {
        newNumbers[j] *= -1; 
      }

      long long newResidue = 0;
      for (auto number : newNumbers) {
              newResidue += number; 
          }
      newResidue = abs(newResidue);

      if (newResidue < currentResidue){
        currentResidue = newResidue;
        currentNumbers = newNumbers;
      }
      else if (exp(-(newResidue - currentResidue) / T(iter)) > distribution(rng)){
        currentNumbers = newNumbers;
        currentResidue = newResidue;
      }
      if (currentResidue < bestResidue)
      {
        bestResidue = currentResidue;
        currentNumbers = newNumbers;
      }
    }

  }
  
  return bestResidue;
}
