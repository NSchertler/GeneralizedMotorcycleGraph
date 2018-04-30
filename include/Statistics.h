#pragma once

#include <vector>

class Statistics
{
public:
	void AddDatum(float datum);

	void Clear();

	std::size_t Size() const;

	void SetGlobalShift(float shift);

	float Average() const;
	float Median();

	//p between 0 and 1
	float Percentile(float p);

	float Variance() const;

	float Min();

	float Max();

private:

	void Sort();

	std::vector<float> data;

	float sum = 0;

	bool isSorted = true;

	float globalShift = 0; //an offset that is added to all data
};