#pragma once

#include "common.h"
#include <map>

//Represents a map from vertices to an arbitrary type that is aware of non-manifold vertices. Non-manifold
//vertices are addressed by their outgoing context halfedge instead of the vertex itself (i.e. they are
//treated as multiple separate manifold vertices).
template <typename T>
class ManifoldnessAwareVertexMap
{
public:
	ManifoldnessAwareVertexMap(const HEMesh& mesh) : mesh(&mesh) { }

	//Removes all entries from the map
	void clear()
	{
		manifoldMap.clear();
		nonManifoldMap.clear();
	}

	//Returns a reference to the value stored for a given manifold vertex.
	//Creates the entry if it does not exist.
	T& AccessOrCreateAtManifoldVertex(HEMesh::VertexHandle v)
	{
		assert(mesh->is_manifold(v));
		return manifoldMap[v];
	}

	//Returns a reference to the value stored for a non-manifold vertex addressed
	//by the given context edge. Creates the entry if it does not exist.
	T& AccessOrCreateAtContextEdge(HEMesh::HalfedgeHandle h)
	{
		assert(mesh->is_boundary(h));
		return nonManifoldMap[h];
	}

	//Returns a reference to the value stored for a non-manifold vertex defined
	//by the target vertex of h. Creates the entry if it does not exist.
	T& AccessOrCreateAtNonManifoldToVertex(HEMesh::HalfedgeHandle h)
	{
		assert(!mesh->is_manifold(mesh->to_vertex_handle(h)));
		//find the context edge
		CirculateBackwardUntil<true>(h, *mesh, [](HEMesh::HalfedgeHandle) { return false; });

		return AccessOrCreateAtContextEdge(h);
	}

	//Returns a reference to the value stored for a vertex addressed
	//by the target vertex of h. Creates the entry if it does not exist.
	T& AccessOrCreateAtToVertex(HEMesh::HalfedgeHandle h)
	{
		auto v = mesh->to_vertex_handle(h);

		if (mesh->is_manifold(v))
			return AccessOrCreateAtManifoldVertex(v);
		else
			return AccessOrCreateAtNonManifoldToVertex(h);
	}

	//Tries to find the entry for the given manifold vertex. Returns if the
	//entry exists. Stores a pointer to the entry in out.
	bool TryAccessAtManifoldVertex(HEMesh::VertexHandle v, T*& out)
	{
		assert(mesh->is_manifold(v));
		auto it = manifoldMap.find(v);
		if (it == manifoldMap.end())
			return false;
		else
		{
			out = &it->second;
			return true;
		}
	}

	//Tries to find the entry for the non-manifold vertex represented by its 
	//context halfedge. Returns if the entry exists. Stores a pointer to the entry in out.
	bool TryAccessAtContextEdge(HEMesh::HalfedgeHandle h, T*& out)
	{
		assert(mesh->is_boundary(h));
		auto it = nonManifoldMap.find(h);
		if (it == nonManifoldMap.end())
			return false;
		else
		{
			out = &it->second;
			return true;
		}
	}

	//Tries to find the entry for the non-manifold vertex represented by one of its 
	//incident halfedges. Returns if the entry exists. Stores a pointer to the entry in out.
	bool TryAccessAtNonManifoldToVertex(HEMesh::HalfedgeHandle h, T*& out)
	{
		assert(!mesh->is_manifold(mesh->to_vertex_handle(h)));
		//find the context edge
		CirculateBackwardUntil<true>(h, *mesh, [](HEMesh::HalfedgeHandle) { return false; });

		return TryAccessAtContextEdge(h, out);
	}

	//Tries to find the entry for the vertex represented by one of its incident 
	//halfedges. Returns if the entry exists. Stores a pointer to the entry in out.
	bool TryAccessAtToVertex(HEMesh::HalfedgeHandle h, T*& out)
	{
		auto v = mesh->to_vertex_handle(h);

		if (mesh->is_manifold(v))
			return TryAccessAtManifoldVertex(v, out);
		else
			return TryAccessAtNonManifoldToVertex(h, out);
	}

	//Tries to find the entry for the given manifold vertex. Returns if the
	//entry exists. Stores a pointer to the entry in out.
	bool TryAccessAtManifoldVertex(HEMesh::VertexHandle v, const T*& out) const
	{
		assert(mesh->is_manifold(v));
		auto it = manifoldMap.find(v);
		if (it == manifoldMap.end())
			return false;
		else
		{
			out = &it->second;
			return true;
		}
	}

	//Tries to find the entry for the non-manifold vertex represented by its 
	//context halfedge. Returns if the entry exists. Stores a pointer to the entry in out.
	bool TryAccessAtContextEdge(HEMesh::HalfedgeHandle h, const T*& out) const
	{
		assert(mesh->is_boundary(h));
		auto it = nonManifoldMap.find(h);
		if (it == nonManifoldMap.end())
			return false;
		else
		{
			out = &it->second;
			return true;
		}
	}

	//Tries to find the entry for the non-manifold vertex represented by one of its 
	//incident halfedges. Returns if the entry exists. Stores a pointer to the entry in out.
	bool TryAccessAtNonManifoldToVertex(HEMesh::HalfedgeHandle h, const T*& out) const
	{
		assert(!mesh->is_manifold(mesh->to_vertex_handle(h)));
		//find the context edge
		CirculateBackwardUntil<true>(h, *mesh, [](HEMesh::HalfedgeHandle) { return false; });

		return TryAccessAtContextEdge(h, out);
	}

	//Tries to find the entry for the vertex represented by one of its incident 
	//halfedges. Returns if the entry exists. Stores a pointer to the entry in out.
	bool TryAccessAtToVertex(HEMesh::HalfedgeHandle h, const T*& out) const
	{
		auto v = mesh->to_vertex_handle(h);

		if (mesh->is_manifold(v))
			return TryAccessAtManifoldVertex(v, out);
		else
			return TryAccessAtNonManifoldToVertex(h, out);
	}

	const T& AccessAtToVertex(HEMesh::HalfedgeHandle h) const
	{
		const T* ptr;
		auto exists = TryAccessAtToVertex(h, ptr);
		if (!exists)
			throw std::runtime_error("The entry for halfedge " + std::to_string(h.idx()) + " does not exist.");
		return *ptr;
	}

	//Returns a reference to the value stored for a given vertex represented by
	//one of its incident halfedges. Throws an exception the entry if it does not exist.
	T& AccessAtToVertex(HEMesh::HalfedgeHandle h)
	{
		T* ptr;
		auto exists = TryAccessAtToVertex(h, ptr);
		if (!exists)
			throw std::runtime_error("The entry for halfedge " + std::to_string(h.idx()) + " does not exist.");
		return *ptr;
	}

	//Erases an etry for a given manifold vertex. Throws an exception if the
	//entry does not exist.
	void EraseAtManifoldVertex(HEMesh::VertexHandle v)
	{
		assert(mesh->is_manifold(v));
		auto it = manifoldMap.find(v);
		if (it == manifoldMap.end())
			throw std::runtime_error("Cannot erase entry for vertex because it does not exist.");
		else
			manifoldMap.erase(it);			
	}

	//Erases an etry for a non-manifold vertex represented by its context edge.
	//Throws an exception if the entry does not exist.
	void EraseAtContextEdge(HEMesh::HalfedgeHandle h)
	{
		assert(mesh->is_boundary(h));
		auto it = nonManifoldMap.find(h);
		if (it == nonManifoldMap.end())
			throw std::runtime_error("Cannot erase entry for context edge because it does not exist.");
		else
			nonManifoldMap.erase(it);
	}

	//Erases an etry for a non-manifold vertex represented by one of its incident halfedges.
	//Throws an exception if the entry does not exist.
	void EraseAtNonManifoldToVertex(HEMesh::HalfedgeHandle h)
	{
		assert(!mesh->is_manifold(mesh->to_vertex_handle(h)));
		//find the context edge
		CirculateBackwardUntil<true>(h, *mesh, [](HEMesh::HalfedgeHandle) { return false; });

		EraseAtContextEdge(h);
	}

	//Erases an etry for a vertex represented by one of its incident halfedges.
	//Throws an exception if the entry does not exist.
	void EraseAtToVertex(HEMesh::HalfedgeHandle h)
	{
		auto v = mesh->to_vertex_handle(h);

		if (mesh->is_manifold(v))
			EraseAtManifoldVertex(v);
		else
			EraseAtNonManifoldToVertex(h);
	}

	//Iterates all entries (manifold and non-manifold) of the ManifoldnessAwareVertexMap.
	class Iterator : public std::iterator<std::forward_iterator_tag, std::pair<HEMesh::HalfedgeHandle, T>>
	{
	public:

		enum State
		{
			Manifold,
			NonManifold
		};

		Iterator()
			: container(nullptr)
		{ }

		Iterator(const ManifoldnessAwareVertexMap<T>& container, State state, 
			typename std::map<HEMesh::VertexHandle, T>::const_iterator manifoldIt, typename std::map<HEMesh::HalfedgeHandle, T>::const_iterator nonManifoldIt)
			: container(&container), state(state), manifoldIt(manifoldIt), nonManifoldIt(nonManifoldIt)
		{ 
			MakeValid();
		}

		Iterator& operator++()
		{
			assert(container != nullptr);
			if (state == Manifold)
			{
				++manifoldIt;
				MakeValid();
			}
			else
			{
				++nonManifoldIt;
			}

			return *this;
		}

		bool operator==(const Iterator& rhs) const
		{
			assert(container != nullptr);
			if (container != rhs.container || state != rhs.state)
				return false;
			if (state == Manifold)
				return manifoldIt == rhs.manifoldIt;
			else
				return nonManifoldIt == rhs.nonManifoldIt;
		}

		bool operator!=(const Iterator& rhs) const { return !(*this == rhs); }

		//Returns a pair of the context outgoing edge of a vertex and the respective stored value.
		typename Iterator::value_type operator*()
		{
			assert(container != nullptr);
			if (state == Manifold)
			{
				auto& entry = *manifoldIt;
				auto contextEdge = container->mesh->halfedge_handle(entry.first);
				return std::make_pair(contextEdge, entry.second);
			}
			else
			{
				auto& entry = *nonManifoldIt;
				return entry;
			}
		}		

	private:
		void MakeValid()
		{
			if (state == Manifold && manifoldIt == container->manifoldMap.end())
			{
				state = NonManifold;
				nonManifoldIt = container->nonManifoldMap.begin();
			}
		}

		const ManifoldnessAwareVertexMap<T>* container;

		State state;

		typename std::map<HEMesh::VertexHandle, T>::const_iterator manifoldIt;
		typename std::map<HEMesh::HalfedgeHandle, T>::const_iterator nonManifoldIt;
	};

	Iterator begin() const { return Iterator(*this, Iterator::State::Manifold, manifoldMap.cbegin(), nonManifoldMap.cbegin()); }
	Iterator end() const { return Iterator(*this, Iterator::State::NonManifold, manifoldMap.cend(), nonManifoldMap.cend()); }

private:
	friend class Iterator;

	std::map<HEMesh::VertexHandle, T> manifoldMap;
	std::map<HEMesh::HalfedgeHandle, T> nonManifoldMap;

	const HEMesh* mesh;
};