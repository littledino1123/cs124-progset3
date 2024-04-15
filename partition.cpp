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

// Define the function to implement the Karmarkar-Karp algorithm.
long long karmarkarKarp(const std::vector<long long> &numbers);

// Define the function for repeated random algorithm.
long long repeatedRandom(const std::vector<long long> &numbers, bool prepartition);

// Define the function for hill climbing algorithm.
long long hillClimbing(const std::vector<long long> &numbers, bool prepartition);

// Define the function for simulated annealing algorithm.
long long simulatedAnnealing(const std::vector<long long> &numbers, bool prepartition);

// Main routine that selects and runs the specified algorithm.
int main(int argc, char *argv[])
{
  if (argc != 4)
  {
    std::cerr << "Usage: " << argv[0] << " flag algorithm inputfile" << std::endl;
    return 1;
  }

  int flag = std::atoi(argv[1]);
  int algorithm = std::atoi(argv[2]);
  std::string inputfile = argv[3];

  // Read the integers from the input file.
  std::ifstream infile(inputfile);
  if (!infile)
  {
    std::cerr << "Error opening file: " << inputfile << std::endl;
    return 1;
  }

  std::vector<long long> numbers;
  long long number;
  while (infile >> number)
  {
    numbers.push_back(number);
  }

  // Check if the input has the correct number of integers.
  if (numbers.size() != 100)
  {
    std::cerr << "Error: Input file should contain 100 integers." << std::endl;
    return 1;
  }

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

long long repeatedRandom(const std::vector<long long> &numbers, bool prepartition)
{

  long long bestResidue = INT_MAX; // You might want to compute this properly if prepartitioning

  for (int iter = 0; iter < 10000; ++iter)
  {
    std::vector<long long> newNumbers = numbers;

    // Generate a random solution
    for (auto &number : newNumbers)
    {
      float w = distribution(rng);
      if (w > 0.5)
      {
        number = -number; // Flip the sign randomly
      }
    }
    long long newResidue = 0;
    for (auto number : newNumbers) {
            newResidue += number; // Add each number to the sum
        }
    bestResidue = std::min(abs(bestResidue), abs(newResidue));
  }

  return bestResidue;
}

long long hillClimbing(const std::vector<long long> &numbers, bool prepartition)
{
  srand(time(NULL));
  long long bestResidue = 0;

  std::vector<long long> currentNumbers = numbers;
  for (auto number : currentNumbers) {
            bestResidue += number; // Add each number to the sum
        }
  for (int iter = 0; iter < 10000; ++iter)
  {
    std::vector<long long> newNumbers = currentNumbers;
    newNumbers[rand() % newNumbers.size()] *= -1; // Flip one element's sign

    long long newResidue = 0;
    for (auto number : newNumbers) {
            newResidue += number; // Add each number to the sum
        }
    newResidue = abs(newResidue);
    if (newResidue < bestResidue)
    {
      bestResidue = newResidue;
      currentNumbers = newNumbers;
    }
  }

  return bestResidue;
}

long long simulatedAnnealing(const std::vector<long long> &numbers, bool prepartition)
{
  srand(time(NULL));
  long long bestResidue = karmarkarKarp(numbers);

  std::vector<long long> currentNumbers = numbers;
  long long currentResidue = bestResidue;

  for (int iter = 0; iter < 10000; ++iter)
  {
    std::vector<long long> newNumbers = currentNumbers;
    newNumbers[rand() % newNumbers.size()] *= -1; // Flip one element's sign

    long long newResidue = karmarkarKarp(newNumbers);
    if (newResidue < currentResidue || exp(-(newResidue - currentResidue) / (double)(iter + 1)) > (rand() / (RAND_MAX + 1.0)))
    {
      currentNumbers = newNumbers;
      currentResidue = newResidue;
    }

    if (currentResidue < bestResidue)
    {
      bestResidue = currentResidue;
    }
  }

  return bestResidue;
}
