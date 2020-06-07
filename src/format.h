//
// Created by xgallom on 4/21/20.
//

#ifndef QM_MONTE_CARLO_SRC_FORMAT_H
#define QM_MONTE_CARLO_SRC_FORMAT_H

#include <string>
#include <ostream>

namespace format {
	extern size_t currentOffset;
	static constexpr size_t OffsetDelta = 4;

	template<typename Member>
	inline void PrintMember(std::ostream &out, const Member &member, const char *name, const std::string &offset)
	{
		out << offset << name << ": " << member << ",\n";
	}

	template<typename Format>
	inline std::ostream &PrintObject(std::ostream &out, const typename Format::Object &object)
	{
		out << Format::Name << "{\n";

		currentOffset += OffsetDelta;
		std::string offset(currentOffset, ' ');

		Format::PrintMembers(out, object, offset);

		currentOffset -= OffsetDelta;
		offset.resize(currentOffset, ' ');

		return out << offset << "}";
	}
}

#endif //QM_MONTE_CARLO_SRC_FORMAT_H
