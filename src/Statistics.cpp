#include "Statistics.h"

#include <iostream>
#include <algorithm>
#include <cmath>

void Statistics::AddDatum(float datum)
{
	if (std::isnan(datum) || std::isinf(datum))
		return;
	if (data.size() > 0 && datum < data.back())
		isSorted = false;
	data.push_back(datum);
	sum += datum;
}

void Statistics::Clear()
{
	sum = 0;
	data.clear();
	isSorted = true;
	globalShift = 0;
}

std::size_t Statistics::Size() const
{
	return data.size();
}

void Statistics::SetGlobalShift(float shift)
{
	globalShift = shift;
}

float Statistics::Average() const
{
	return sum / Size() + globalShift;
}
float Statistics::Median()
{
	return Percentile(0.5f);
}

//p between 0 and 1
float Statistics::Percentile(float p)
{
	if (Size() == 0)
		return std::numeric_limits<float>::quiet_NaN();
	Sort();
	float index = p * (Size() - 1);

	auto lower = std::floor(index);
	auto lowerWeight = 1 + lower - index;

	float percentile = lowerWeight * data[(int)lower];
	if (lowerWeight != 1)
		percentile += (1 - lowerWeight) * data[(int)lower + 1];

	return percentile + globalShift;
}

float Statistics::Variance() const
{
	auto avg = Average();

	float variance = 0;
	for (auto r = data.begin(); r != data.end(); ++r)
	{
		float d = (*r + globalShift) - avg;
		variance += d * d;
	}
	return variance / (Size() - 1);
}

float Statistics::Min()
{
	Sort();
	return data[0];
}

float Statistics::Max()
{
	Sort();
	return data.back();
}

void Statistics::Sort()
{
	if (isSorted)
		return;
	std::sort(data.begin(), data.end());
	isSorted = true;
}
