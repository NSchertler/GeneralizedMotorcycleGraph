#include "PolygonTriangulation.h"

#include <functional>

float SaneArg(float v)
{
	return std::max(-1.0f, std::min(v, 1.0f));
}

//Maximize minimal interior angle
float TriangleCost(const Eigen::Vector3f& v1, const Eigen::Vector3f& v2, const Eigen::Vector3f& v3)
{
	float lSq1 = (v2 - v1).squaredNorm();
	float lSq2 = (v3 - v2).squaredNorm();
	float lSq3 = (v1 - v3).squaredNorm();

	float l1 = sqrt(lSq1);
	float l2 = sqrt(lSq2);
	float l3 = sqrt(lSq3);

	float angle1 = std::acos(SaneArg((lSq2 + lSq3 - lSq1) / (2 * l2 * l3)));
	float angle2 = std::acos(SaneArg((lSq3 + lSq1 - lSq2) / (2 * l3 * l1)));
	float angle3 = std::acos(SaneArg((lSq1 + lSq2 - lSq3) / (2 * l1 * l2)));

	float minAngle;
	if (angle1 < angle2 && angle1 < angle3)
		minAngle = angle1;
	else if (angle2 < angle3)
		minAngle = angle2;
	else
		minAngle = angle3;
	return -minAngle;
}

std::vector<std::array<int, 3>> TriangulatePolygon(const std::vector<Eigen::Vector3f>& polygon)
{
	if (polygon.size() < 3)
		return std::vector<std::array<int, 3>>();

	//variable names are chosen according to the paper
	struct DPEntry
	{
		float weight; //W
		int minimumIndex; //O
	};
	std::vector<std::vector<DPEntry>> W_ji(polygon.size());
	for (int j = 1; j < polygon.size(); ++j)
		W_ji[j].resize(j);

	//algorithm step 1
	for (int i = 0; i <= polygon.size() - 2; ++i)
		W_ji[i + 1][i].weight = 0;
	for (int i = 0; i <= polygon.size() - 3; ++i)
		W_ji[i + 2][i].weight = TriangleCost(polygon[i], polygon[i + 1], polygon[i + 2]);

	//algorithm step 2, 3
	for (int j = 3; j < polygon.size(); ++j)
	{
		for (int i = 0; i <= polygon.size() - j - 1; ++i)
		{
			int k = i + j;
			auto& entry = W_ji[k][i];

			float minCost = std::numeric_limits<float>::infinity();
			int minIndex;
			for (int m = i + 1; m < k; ++m)
			{
				float cost = W_ji[m][i].weight + W_ji[k][m].weight + TriangleCost(polygon[i], polygon[m], polygon[k]);
				if (cost < minCost)
				{
					minCost = cost;
					minIndex = m;
				}
			}

			entry.minimumIndex = minIndex;
			entry.weight = minCost;
		}
	}

	std::vector<std::array<int, 3>> resultFaces;
	std::function<void(int, int)> trace = [&](int i, int k)
	{
		if (i + 2 == k)
			resultFaces.push_back({ i, i + 1, k });
		else
		{
			int o = W_ji[k][i].minimumIndex;
			if (o != i + 1)
				trace(i, o);
			resultFaces.push_back({ i, o, k });
			if (o != k - 1)
				trace(o, k);
		}
	};

	//algorithm step 4
	trace(0, polygon.size() - 1);

	return resultFaces;
}