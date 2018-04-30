/*
	@author Wenzel Jakob
	@author Nico Schertler
*/

#include "meshio.h"
#include <nsessentials/data/FileHelper.h>
#include <fstream>
#include <memory>

#include <nsessentials/util/Timer.h>

extern "C" {
    #include "rply.h"
}

void load_ply(const std::string &filename, Matrix3Xf &V,
	FaceList& F)
{
	//Matrix3Xf &N, Matrix3Xus &C;

	auto message_cb = [](p_ply ply, const char *msg) { std::cerr << "rply: " << msg << std::endl; };

	nse::util::Timer<> timer;
	std::cout << "Loading \"" << filename << "\" .. ";
	std::cout.flush();

	p_ply ply = ply_open(filename.c_str(), message_cb, 0, nullptr);
	if (!ply)
		throw std::runtime_error("Unable to open PLY file \"" + filename + "\"!");

	if (!ply_read_header(ply)) {
		ply_close(ply);
		throw std::runtime_error("Unable to open PLY header of \"" + filename + "\"!");
	}

	const float gamma = 2.2f;

	p_ply_element element = nullptr;
	uint32_t vertexCount = 0, faceCount = 0;

	p_ply_property prop = nullptr;
	const char* pname;
	e_ply_type type, length, value;

	bool hasUv = false;
	bool hasNormals = false;
	bool hasColor = false;

	/* Inspect the structure of the PLY file, load number of faces if avaliable */
	while ((element = ply_get_next_element(ply, element)) != nullptr) {
		const char *name;
		long nInstances;

		ply_get_element_info(element, &name, &nInstances);
		if (!strcmp(name, "vertex"))
		{
			vertexCount = (uint32_t)nInstances;
			while ((prop = ply_get_next_property(element, prop)) != nullptr)
			{
				ply_get_property_info(prop, &pname, &type, &length, &value);
				if (!strcmp(pname, "texture_u"))
					hasUv = true;
				else if (!strcmp(pname, "nx"))
					hasNormals = true;
				else if (!strcmp(pname, "red"))
					hasColor = true;
			}
		}
		else if (!strcmp(name, "face"))
			faceCount = (uint32_t)nInstances;
	}

	if (vertexCount == 0 || faceCount == 0)
		throw std::runtime_error("PLY file \"" + filename + "\" is invalid! No face/vertex/elements found!");

	F.resize(faceCount);
	V.resize(3, vertexCount);
	//N.resize(3, hasNormals ? vertexCount : 0);
	//C.resize(3, hasColor ? vertexCount : 0);
	//Eigen::Matrix<float, 2, Eigen::Dynamic> uv(2, hasUv ? vertexCount : 0);


	struct VertexCallbackData {
		Matrix3Xf &V;
		VertexCallbackData(Matrix3Xf &V)
			: V(V) { }
	};

	//struct VertexNormalCallbackData {
	//	Matrix3Xf &N;
	//	VertexNormalCallbackData(Matrix3Xf &_N)
	//		: N(_N) { }
	//};
	//
	//struct VertexUVCallbackData {
	//	Eigen::Matrix<float, 2, Eigen::Dynamic> &uv;
	//	VertexUVCallbackData(Eigen::Matrix<float, 2, Eigen::Dynamic> &_uv)
	//		: uv(_uv) { }
	//};
	//
	//struct VertexColorCallbackData {
	//	Eigen::Matrix<unsigned short, 3, Eigen::Dynamic> &c;
	//	VertexColorCallbackData(Eigen::Matrix<unsigned short, 3, Eigen::Dynamic> &c)
	//		: c(c) { }
	//};

	auto rply_vertex_cb = [](p_ply_argument argument) -> int {
		VertexCallbackData *data; long index, coord;
		ply_get_argument_user_data(argument, (void **)&data, &coord);
		ply_get_argument_element(argument, nullptr, &index);
		data->V(coord, index) = (float)ply_get_argument_value(argument);
		return 1;
	};

	//auto rply_vertex_normal_cb = [](p_ply_argument argument) -> int {
	//	VertexNormalCallbackData *data; long index, coord;
	//	ply_get_argument_user_data(argument, (void **)&data, &coord);
	//	ply_get_argument_element(argument, nullptr, &index);
	//	data->N(coord, index) = (float)ply_get_argument_value(argument);
	//	return 1;
	//};
	//
	//auto rply_vertex_uv_cb = [](p_ply_argument argument) -> int
	//{
	//	VertexUVCallbackData *data; long index, coord;
	//	ply_get_argument_user_data(argument, (void **)&data, &coord);
	//	ply_get_argument_element(argument, nullptr, &index);
	//	data->uv(coord, index) = (float)ply_get_argument_value(argument);
	//	return 1;
	//};
	//
	//auto rply_vertex_color_cb = [](p_ply_argument argument) -> int
	//{
	//	VertexColorCallbackData *data; long index, coord;
	//	ply_get_argument_user_data(argument, (void **)&data, &coord);
	//	ply_get_argument_element(argument, nullptr, &index);
	//	auto colorInPly = ply_get_argument_value(argument);
	//	data->c(coord, index) = (unsigned short)(std::pow(colorInPly / 255.0, 2.2) * 65535);
	//	return 1;
	//};

	auto rply_index_cb = [](p_ply_argument argument) -> int {
		FaceList *data;
		long length, value_index, index;
		ply_get_argument_property(argument, nullptr, &length, &value_index);

		ply_get_argument_user_data(argument, (void **)&data, nullptr);
		ply_get_argument_element(argument, nullptr, &index);

		auto& faceEntry = data->at(index);
		if (faceEntry.size() == 0)
			faceEntry.resize(length);

		if (value_index >= 0)
			faceEntry.at(value_index) = (uint32_t)ply_get_argument_value(argument);

		return 1;
	};

	VertexCallbackData vcbData(V);
	//VertexNormalCallbackData vncbData(N);
	//VertexUVCallbackData vuvcbData(uv);
	//VertexColorCallbackData vcData(C);

	ply_set_read_cb(ply, "vertex", "x", rply_vertex_cb, &vcbData, 0);
	ply_set_read_cb(ply, "vertex", "y", rply_vertex_cb, &vcbData, 1);
	ply_set_read_cb(ply, "vertex", "z", rply_vertex_cb, &vcbData, 2);

	//ply_set_read_cb(ply, "vertex", "nx", rply_vertex_normal_cb, &vncbData, 0);
	//ply_set_read_cb(ply, "vertex", "ny", rply_vertex_normal_cb, &vncbData, 1);
	//ply_set_read_cb(ply, "vertex", "nz", rply_vertex_normal_cb, &vncbData, 2);

	if (faceCount > 0)
	{
		if (!ply_set_read_cb(ply, "face", "vertex_index", rply_index_cb, &F, 0) && !ply_set_read_cb(ply, "face", "vertex_indices", rply_index_cb, &F, 0))
		{
			ply_close(ply);
			throw std::runtime_error("PLY file \"" + filename + "\" does not contain vertex indices!");
		}
	}

	//ply_set_read_cb(ply, "vertex", "texture_u", rply_vertex_uv_cb, &vuvcbData, 0);
	//ply_set_read_cb(ply, "vertex", "texture_v", rply_vertex_uv_cb, &vuvcbData, 1);
	//
	//ply_set_read_cb(ply, "vertex", "red", rply_vertex_color_cb, &vcData, 0);
	//ply_set_read_cb(ply, "vertex", "green", rply_vertex_color_cb, &vcData, 1);
	//ply_set_read_cb(ply, "vertex", "blue", rply_vertex_color_cb, &vcData, 2);

	if (!ply_read(ply)) {
		ply_close(ply);
		throw std::runtime_error("Error while loading PLY data from \"" + filename + "\"!");
	}

	ply_close(ply);

	std::cout << "done. (V=" << vertexCount;
	if (faceCount > 0)
		std::cout << ", F=" << faceCount;
	std::cout << ", took " << nse::util::timeString(timer.value()) << ")" << std::endl;
}

inline uint32_t str_to_uint32_t(const std::string &str) {
	char *end_ptr = nullptr;
	uint32_t result = (uint32_t)strtoul(str.c_str(), &end_ptr, 10);
	if (*end_ptr != '\0')
		throw std::runtime_error("Could not parse unsigned integer \"" + str + "\"");
	return result;
}

inline uint32_t str_to_int32_t(const std::string &str) {
	char *end_ptr = nullptr;
	int32_t result = (int32_t)strtol(str.c_str(), &end_ptr, 10);
	if (*end_ptr != '\0')
		throw std::runtime_error("Could not parse signed integer \"" + str + "\"");
	return result;
}

inline float str_to_float(const std::string &str) {
	char *end_ptr = nullptr;
	float result = (float)strtod(str.c_str(), &end_ptr);
	if (*end_ptr != '\0')
		throw std::runtime_error("Could not parse floating point value \"" + str + "\"");
	return result;
}

inline std::vector<std::string> &str_tokenize(const std::string &s, char delim, std::vector<std::string> &elems, bool include_empty = false) {
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim))
		if (!item.empty() || include_empty)
			elems.push_back(item);
	return elems;
}

inline std::vector<std::string> str_tokenize(const std::string &s, char delim, bool include_empty) {
	std::vector<std::string> elems;
	str_tokenize(s, delim, elems, include_empty);
	return elems;
}

void load_obj(const std::string &filename, Matrix3Xf &V, FaceList& F, Matrix2Xf &T, FaceList& FT)
{
	std::ifstream is(filename);
	if (is.fail())
		throw std::runtime_error("Unable to open OBJ file \"" + filename + "\"!");
	std::cout << "Loading \"" << filename << "\" .. ";
	std::cout.flush();
	nse::util::Timer<> timer;


	size_t nPositions = 0;
	V.resize(3, 32);
	size_t nTex = 0;
	T.resize(2, 32);

	std::string line_str;
	while (std::getline(is, line_str)) {
		std::istringstream line(line_str);

		std::string prefix;
		line >> prefix;

		if (prefix == "v") {
			if (nPositions + 1 > V.cols())
				V.conservativeResize(3, 2 * V.cols());
			line >> V.coeffRef(0, nPositions) >> V.coeffRef(1, nPositions) >> V.coeffRef(2, nPositions);
			++nPositions;
		}
		if (prefix == "vt") {
			if (nTex + 1 > T.cols())
				T.conservativeResize(2, 2 * T.cols());
			line >> T.coeffRef(0, nTex) >> T.coeffRef(1, nTex) >> T.coeffRef(2, nTex);
			++nTex;
		}
		else if (prefix == "f") {
			std::string v;
			F.emplace_back();
			if (nTex > 0)
				FT.emplace_back();
			while (line >> v)
			{
				std::vector<std::string> tokens = str_tokenize(v, '/', true);

				if (tokens.size() < 1 || tokens.size() > 3)
					throw std::runtime_error("Invalid vertex data: \"" + v + "\"");

				auto vIdx = str_to_uint32_t(tokens[0]);
				F.back().push_back(vIdx - 1);

				if (tokens.size() >= 2 && !tokens[1].empty())
				{
					auto tIdx = str_to_uint32_t(tokens[1]);
					FT.back().push_back(tIdx - 1);
				}
			}			
		}
	}	

	V.conservativeResize(3, nPositions);	
	T.conservativeResize(2, nTex);

	std::cout << "done. (V=" << V.cols() << ", F=" << F.size() << ", took "
		<< nse::util::timeString(timer.value()) << ")" << std::endl;
}

void load_obj(const std::string &filename, Matrix3Xf &V, FaceList& F)
{
	Matrix2Xf T;
	FaceList FT;
	load_obj(filename, V, F, T, FT);
}

void LoadMesh(const std::string &filename, Matrix3Xf& V, FaceList& F)
{
	std::string ext = nse::data::extension(filename);

	std::string cachePath = filename + ".cache";

	/*if (nse::data::file_exists(cachePath))
	{
		//load from cache

		FILE* f = fopen(cachePath.c_str(), "rb");
		size_t n;
		fread(&n, sizeof(size_t), 1, f);
		V.resize(3, n);
		fread(V.data(), sizeof(float), V.size(), f);

		fread(&n, sizeof(size_t), 1, f);
		N.resize(3, n);
		fread(N.data(), sizeof(float), N.size(), f);

		fread(&n, sizeof(size_t), 1, f);
		C.resize(3, n);
		fread(C.data(), sizeof(unsigned short), C.size(), f);

		fread(&n, sizeof(size_t), 1, f);
		F.resize(3, n);
		fread(F.data(), sizeof(uint32_t), F.size(), f);
		fclose(f);

		//scan = new Scan(V, N, C, F, nse::data::filename_without_extension_and_directory(filename), Eigen::Affine3f(transform));
	}
	else*/
	{
		if (ext == ".ply")
			load_ply(filename, V, F);
		else if (ext == ".obj")
			load_obj(filename, V, F);
		else
		{
			throw std::runtime_error("LoadMesh: Unknown file extension.");
		}

		//scan = new Scan(V, N, C, F, nse::data::filename_without_extension_and_directory(filename), Eigen::Affine3f(transform));

		/*FILE* f = fopen(cachePath.c_str(), "wb");
		size_t n = V.cols();
		fwrite(&n, sizeof(size_t), 1, f);
		fwrite(V.data(), sizeof(float), V.size(), f);
		n = N.cols();
		fwrite(&n, sizeof(size_t), 1, f);
		fwrite(N.data(), sizeof(float), N.size(), f);
		n = C.cols();
		fwrite(&n, sizeof(size_t), 1, f);
		fwrite(C.data(), sizeof(unsigned short), C.size(), f);
		n = F.cols();
		fwrite(&n, sizeof(size_t), 1, f);
		fwrite(F.data(), sizeof(uint32_t), F.size(), f);
		fclose(f);*/
	}

	//dataSink.AddScan(scan);		
}