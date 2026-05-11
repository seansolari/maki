
#include "maki/build/graph/construct_cdbg.hpp"
#include "maki/build/graph/construct_common.hpp"
#include "maki/core/graph/cdbg.hpp"
#include "maki/core/graph/wdbg.hpp"
#include "maki/core/seq/io.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(_maki, m) {

  // graph interface

  py::class_<ColouredGraphFiles>(m, "cdbg_files").def(py::init<fs::path>());

  py::class_<ColouredGraph>(m, "cdbg")
      .def(py::init<>())
      .def_static("disk_load", &ColouredGraph::FromDisk,
                  "Load a coloured graph from disk.");

  py::class_<WeightedGraphFiles>(m, "wdbg_files").def(py::init<fs::path>());

  py::class_<WeightedGraph>(m, "wdbg")
      .def(py::init<>())
      .def_static("disk_load", &WeightedGraph::FromDisk,
                  "Load a coloured graph from disk.");

  // build helpers

  py::enum_<InputFileType>(m, "FileType")
      .value("FNA", InputFileType::FastaFileType)
      .value("GFF3", InputFileType::Gff3FileType)
      .value("FQ", InputFileType::FastQFileType)
      .export_values();

  py::class_<GenomeManifest>(m, "genome_manifest")
      .def(py::init<>())
      .def_readonly("files", &GenomeManifest::files)
      .def_readonly("type", &GenomeManifest::type);

  m.def("read_manifest", &readFilePaths);

  py::class_<dbg::BuildOptions>(m, "build_opts")
      .def(py::init<>())
      .def_readwrite("kmer_size", &dbg::BuildOptions::kmer_size)
      .def_readwrite("suffix_size", &dbg::BuildOptions::suffix_size)
      .def_readwrite("out", &dbg::BuildOptions::out)
      .def_readwrite("pool_size", &dbg::BuildOptions::pool_size)
      .def_readwrite("reserve_per_chunk", &dbg::BuildOptions::reserve_per_chunk)
      .def_readwrite("threads", &dbg::BuildOptions::threads);

  // build API

  m.def("construct_cdbg",
        py::overload_cast<const GenomeManifest &, dbg::BuildOptions>(
            &cdbg::construct));

}
