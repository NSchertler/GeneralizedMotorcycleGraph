#include "Singularity.h"

const std::vector<HEMesh::HalfedgeHandle>& MetaSingularity::GetEmanatingMotorcycle(int orientation) const
{
	return emanatingMotorcycles[orientation].first;
}

std::vector<HEMesh::HalfedgeHandle>& MetaSingularity::GetEmanatingMotorcycle(int orientation)
{
	if (orientation >= emanatingMotorcycles.size())
		emanatingMotorcycles.resize(orientation + 1);
	return emanatingMotorcycles[orientation].first;
}

void MetaSingularity::SetEmanatingMotorcycleActive(int orientation, bool active)
{
	emanatingMotorcycles[orientation].second = active;
}

bool MetaSingularity::GetEmanatingMotorcycleActive(int orientation) const
{
	return emanatingMotorcycles[orientation].second;
}

size_t MetaSingularity::FindOutgoingMotorcycle(HEMesh::HalfedgeHandle h) const
{
	for (int i = 0; i < emanatingMotorcycles.size(); ++i)
		if (!emanatingMotorcycles[i].first.empty() && emanatingMotorcycles[i].first.front() == h)
			return i;
	return (size_t)-1;
}

size_t MetaSingularity::Degree() const
{
	return emanatingMotorcycles.size();
}

void MetaSingularity::AddMotorcycle(HEMesh::HalfedgeHandle direction, bool active)
{
	emanatingMotorcycles.emplace_back();
	emanatingMotorcycles.back().first.push_back(direction);
	emanatingMotorcycles.back().second = active;
}

void MetaSingularity::AddEmptyMotorcycle()
{
	emanatingMotorcycles.emplace_back();
}

void MetaSingularity::AddMotorcycleAtFront(HEMesh::HalfedgeHandle direction, bool active)
{
	emanatingMotorcycles.emplace(emanatingMotorcycles.begin());
	emanatingMotorcycles.front().first.push_back(direction);
	emanatingMotorcycles.front().second = active;
}
