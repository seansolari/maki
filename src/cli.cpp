#include <glog/logging.h>
#include <oneapi/tbb/global_control.h>
#include <oneapi/tbb/enumerable_thread_specific.h>
#include <filesystem>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <syncstream>
#include <maki/cli.hpp>
#include <maki/graph.hpp>
#include <maki/fasta.hpp>
#include <maki/classify.hpp>
#include <maki/dist.hpp>
#include <maki/utils.hpp>
#include <maki/timing.hpp>
#include <maki/maki.h>

using namespace std;
namespace fs = filesystem;

int BaseCommand::run(const char *execName, seqan3::argument_parser &parser)
{
  registerOptions(parser);
  try
  {
    parser.parse();
  }
  catch (seqan3::argument_parser_error &e)
  {
    LOG(ERROR) << "Argument parsing failed: " << e.what() << ". Please check your command syntax.";
    return -1;
  }
  LoggingUtils::setupLogging(logFile, execName);
  ThreadUtils::configureThreads(threads);
  try
  {
    return execute();
  }
  catch (const exception &e)
  {
    LOG(ERROR) << "Execution failed: " << e.what();
    return 1;
  }
}

// ============================
// Argument registration
// ============================
void IndexCommand::registerOptions(seqan3::argument_parser &parser)
{
  parser.info.short_description = "index a collection of genome sequences";
  parser.add_option(threads, 't', "threads", "Number of threads to use for indexing");
  parser.add_option(resourceLimit, 'R', "ram-gb", "RAM limit");
  parser.add_flag(useBulk, 'b', "bulk", "Force bulk indexing strategy");
  parser.add_option(suffixSize, 's', "suffix-size", "Force indexing with this suffix size");
  parser.add_option(queryFile, 'q', "qf", "Path to query file containing genomes to index");
  parser.add_option(outputFolder, 'o', "out", "Path to save graph");
  parser.add_option(kmerSize, 'k', "kmer", "K-mer size");
  parser.add_option(logFile, 'l', "log", "File to write logs");
  parser.add_flag(overrideDatabase, 'F', "force", "Force overwrite of database");
}

void MakeFilterCommand::registerOptions(seqan3::argument_parser &parser)
{
  parser.info.short_description = "create a filter database (input is FASTA read from STDIN)";
  parser.add_option(threads, 't', "threads", "Number of threads to use for indexing");
  parser.add_option(resourceLimit, 'R', "ram-gb", "RAM limit");
  parser.add_flag(useBulk, 'b', "bulk", "Force bulk indexing strategy");
  parser.add_option(suffixSize, 's', "suffix-size", "Force indexing with this suffix size");
  parser.add_option(outputFolder, 'o', "out", "Path to save filter graph");
  parser.add_option(kmerSize, 'k', "kmer", "K-mer size");
  parser.add_option(logFile, 'l', "log", "File to write logs");
  parser.add_flag(overrideDatabase, 'F', "force", "Force overwrite of database");
}

void DbStatsCommand::registerOptions(seqan3::argument_parser &parser)
{
  parser.info.short_description = "Print basic index statistics";
  parser.add_option(databasePath, 'i', "db", "Path to genome index");
  parser.add_flag(loadBig, '\0', "long", "Load full database for more stats");
}

void ClassifyCommand::registerOptions(seqan3::argument_parser &parser)
{
  parser.info.short_description = "analyse collection of short reads";
  parser.add_option(threads, 't', "threads", "Number of threads to use for indexing");
  parser.add_option(databasePath, 'i', "db", "Path to genome index");
  parser.add_option(queryPath, 'q', "query", "Path to query database");
  parser.add_option(outputFile, 'o', "out", "Output file");
  parser.add_option(minReadLength, 'L', "min-length", "Minimum read length");
  parser.add_option(resourceLimit, 'R', "ram-gb", "RAM limit");
  parser.add_flag(useBulk, 'b', "bulk", "Force bulk indexing strategy");
  parser.add_option(suffixSize, 's', "suffix-size", "Force indexing with this suffix size");
  parser.add_option(logFile, 'l', "log", "File to write logs");
}

void DistanceCommand::registerOptions(seqan3::argument_parser &parser)
{
  parser.info.short_description = "compute shared k-mers between genomes using a database";
  parser.add_option(databasePath, 'i', "db", "Path to genome index");
  parser.add_option(outputFile, 'o', "out", "File to save results");
  parser.add_option(threads, 't', "threads", "Number of threads to use for distance calculation");
  parser.add_option(logFile, 'l', "log", "File to write logs");
  parser.add_option(tempDir, 'z', "temp-dir", "Folder to use for temporary data");
}

void CountCommand::registerOptions(seqan3::argument_parser &parser)
{
  parser.info.short_description = "count unique k-mers";
  parser.add_option(databasePath, 'i', "db", "Path to genome index");
  parser.add_option(outputFile, 'o', "out", "File to save results");
  parser.add_option(threads, 't', "threads", "Number of threads to use for indexing");
  parser.add_option(logFile, 'l', "log", "File to write logs");
  parser.add_flag(useWhole, 'W', "whole", "Calculate distance between whole genomes, not features");
}

void CoverageCommand::registerOptions(seqan3::argument_parser &parser)
{
  parser.info.short_description = "calculate coverage";
  parser.add_option(inputFile, 'i', "input", "Input TSV keys and counts");
  parser.add_option(countFile, 'c', "count", "Count output");
  parser.add_option(outputFile, 'o', "out", "File to save results");
  parser.add_flag(skipHeader, '\0', "skip-header", "Skip header line");
  parser.add_option(valueColumn, 'X', "value-column", "Column to take k-mer count of colours from");
}

vector<size_t> SubcountCommand::parseSubsamples() const
{
  vector<size_t> result;
  string numstr;
  stringstream is(countString);
  while (getline(is, numstr, ','))
  {
    try
    {
      result.push_back(stoi(numstr));
    }
    catch (const invalid_argument &e)
    {
      LOG(ERROR) << "Invalid argument: " << numstr << " - " << e.what();
    }
    catch (const out_of_range &e)
    {
      LOG(ERROR) << "Out of range: " << numstr << " - " << e.what();
    }
  }
  return result;
}

vector<size_t> SubcountCommand::extractSubsamples(uint64_t maxId) const
{
  auto values = parseSubsamples();
  size_t originalSize = values.size();
  values.erase(remove_if(values.begin(), values.end(), [&](size_t v)
                         { return v > maxId; }),
               values.end());
  if (values.size() < originalSize)
    LOG(WARNING) << "removed " << originalSize - values.size() << " sample sizes above the maximum value of " << maxId;
  return values;
}

void SubcountCommand::registerOptions(seqan3::argument_parser &parser)
{
  parser.info.short_description = "count unique k-mers in subsets of total";
  parser.add_option(databasePath, 'i', "db", "Path to genome index");
  parser.add_option(outputFile, 'o', "out", "File to save results");
  parser.add_option(countString, 'X', "subsamples", "Comma-separated number of genomes to subsample");
  parser.add_option(threads, 't', "threads", "Number of threads to use for indexing");
  parser.add_option(reps, 'r', "repeats", "Number of repeats per subsample size");
  parser.add_option(logFile, 'l', "log", "File to write logs");
}

void CompressCommand::registerOptions(seqan3::argument_parser &parser)
{
  parser.info.short_description = "count oc-occurring k-mers";
  parser.add_option(databasePath, 'i', "db", "Path to genome index");
  parser.add_option(threads, 't', "threads", "Number of threads to use for counting");
  parser.add_option(outputFile, 'o', "out", "File to save results");
  parser.add_option(logFile, 'l', "log", "File to write logs");
}

void JaccardCommand::registerOptions(seqan3::argument_parser &parser)
{
  parser.info.short_description = "Calculate Jaccard index, p-values and metadata";
  parser.add_option(databasePath, 'i', "db", "Path to genome index");
  parser.add_option(nwkFile, 't', "tree", "Taxonomy of genomes");
  parser.add_option(countFile, 'c', "count", "Count output");
  parser.add_option(distFile, 'd', "dist", "Pairwise distance output");
  parser.add_option(take, 'r', "take", "Number of random records to output per rank");
  parser.add_option(alpha, 'a', "alpha", "Significance level for r1 CI");
  parser.add_option(outputFile, 'o', "out", "File to save results");
}

void LcaCommand::registerOptions(seqan3::argument_parser &parser)
{
  parser.info.short_description = "calculate conservation of k-mers";
  parser.add_option(databasePath, 'i', "db", "Path to genome index");
  parser.add_option(countFile, 'c', "count", "Cooc output");
  parser.add_option(clusterFile, 'g', "cluster", "Cluster file contains 2 tab-separated columns: [representative] [member]");
  parser.add_option(nwkFile, 't', "tree", "Taxonomy of genomes");
  parser.add_option(threads, 'p', "threads", "Number of threads to use for counting");
  parser.add_option(outputFile, 'o', "out", "File to save results");
  parser.add_flag(printColours, 'P', "print-colours", "Print colours (increases file size)");
  parser.add_option(minFamDiv, '\0', "min-div", "Minimum gene family diversity");
  parser.add_option(maxFamDiv, '\0', "max-div", "Maximum gene family diversity");
  parser.add_option(minPhyCons, '\0', "min-cons", "Minimum phylogenetic conservation");
  parser.add_option(maxPhyCons, '\0', "max-cons", "Maximum phylogenetic conservation");
}

void GlcaCommand::registerOptions(seqan3::argument_parser &parser)
{
  parser.info.short_description = "calculate conservation of genes";
  parser.add_option(databasePath, 'i', "db", "Path to genome index");
  parser.add_option(clusterFile, 'g', "cluster", "Cluster file contains 2 tab-separated columns: [representative] [member]");
  parser.add_option(nwkFile, 't', "tree", "Taxonomy of genomes");
  parser.add_option(outputFile, 'o', "out", "File to save results");
}

void SamdistCommand::registerOptions(seqan3::argument_parser &parser)
{
  parser.info.short_description = "sample a feature distance among every pair of sequences";
  parser.add_option(databasePath, 'i', "db", "Path to genome index");
  parser.add_option(countFile, 'c', "count", "Count output");
  parser.add_option(nwkFile, 't', "tree", "Taxonomy of genomes");
  parser.add_option(threads, 'p', "threads", "Number of threads to use for counting");
  parser.add_option(outputFile, 'o', "out", "File to save results");
}

void QueryCommand::registerOptions(seqan3::argument_parser &parser)
{
  parser.info.short_description = "retrieve k-mer counts for specific genomes";
  parser.add_option(queryFile, 'q', "qf", "Path to query file containing genome names to search for");
  parser.add_option(databasePath, 'i', "db", "Path to genome index");
  parser.add_option(countFile, 'c', "count", "Count output");
  parser.add_option(nwkFile, 't', "tree", "Taxonomy of genomes");
  parser.add_option(threads, 'p', "threads", "Number of threads to use for counting");
  parser.add_option(outputFile, 'o', "out", "File to save results");
}

void KmersCommand::registerOptions(seqan3::argument_parser &parser)
{
  parser.info.short_description = "retrieve k-mers with specific colour profiles";
  parser.add_option(coloursFile, 'c', "colours", "Path to query file containing colour profiles");
  parser.add_option(databasePath, 'i', "db", "Path to genome index");
  parser.add_option(threads, 'p', "threads", "Number of threads to use for counting");
  parser.add_option(outputFile, 'o', "out", "File to save results");
  parser.add_flag(skipHeader, '\0', "skip-header", "Skip header line");
}

// ============================
// IndexCommand
// ============================
int IndexCommand::execute()
{
  if (!FilesystemUtils::checkPathExists(queryFile, "Query file"))
    return -1;
  if (!ValidationUtils::validateIndexingFlags(useBulk, suffixSize, resourceLimit))
    return -1;
  if (!FilesystemUtils::ensureDirectoryExists(outputFolder, overrideDatabase))
    return -1;
  LOG(INFO) << "Running index on query file=" << queryFile << " using threads=" << threads;
  try
  {
    GraphBuilder::IndexParams params{queryFile, outputFolder, static_cast<uint8_t>(kmerSize - 1u), threads, resourceLimit, useBulk, suffixSize};
    GraphBuilder::buildIndex(params);
  }
  catch (const exception &e)
  {
    LOG(ERROR) << "Indexing failed: " << e.what();
    return 1;
  }
  return 0;
}

// ============================
// MakeFilterCommand
// ============================
int MakeFilterCommand::execute()
{
  if (!FilesystemUtils::ensureDirectoryExists(outputFolder, overrideDatabase))
    return -1;
  LOG(INFO) << "Creating filter database using threads=" << threads;
  try
  {
    GraphBuilder::IndexParams params{"stdin", outputFolder, static_cast<uint8_t>(kmerSize - 1u), threads, resourceLimit, useBulk, suffixSize};
    GraphBuilder::buildFilter(params);
  }
  catch (const exception &e)
  {
    LOG(ERROR) << "Filter creation failed: " << e.what();
    return 1;
  }
  return 0;
}

// ============================
// DbStatsCommand
// ============================
int DbStatsCommand::execute()
{
  if (!FilesystemUtils::checkPathExists(databasePath, "Database"))
    return -1;
  LOG(INFO) << "Inspecting database at " << databasePath;
  try
  {
    graph::TraversableGraph db;
    if (loadBig)
    {
      graph::loadFromDisk(databasePath, db);
      cout << "k-mer size: " << (size_t)db.k
           << "\nmax colour: " << db.getMaxColour()
           << "\nnum genomes: " << db.getNumGenomes()
           << "\ncolour bit width: " << (size_t)db.getColourBitWidth()
           << "\nbase bit width: " << (size_t)db.getBaseColourBitWidth()
           << "\nnum edges: " << db.edge_count()
           << '\n';
    }
    else
    {
      graph::StaticGraphDiskFileConfig config(databasePath);
      graph::loadSmallBuffers(config.meta, db);
      cout << "k-mer size: " << (size_t)db.k
           << "\nmax colour: " << db.getMaxColour()
           << "\nnum genomes: " << db.getNumGenomes()
           << "\ncolour bit width: " << (size_t)db.getColourBitWidth()
           << "\nbase bit width: " << (size_t)db.getBaseColourBitWidth()
           << '\n';
    }
  }
  catch (const exception &e)
  {
    LOG(ERROR) << "Failed to load database: " << e.what();
    return 1;
  }
  return 0;
}

// ============================
// ClassifyCommand
// ============================
int ClassifyCommand::execute()
{
  if (!FilesystemUtils::checkPathExists(databasePath, "Database"))
    return -1;
  if (!FilesystemUtils::checkPathExists(queryPath, "Query file"))
    return -1;
  try
  {
    graph::TraversableGraph db;
    graph::loadFromDisk(databasePath, db);
    LOG(INFO) << "loading filter database from " << queryPath;
    graph::TraversableGraph filter;
    graph::loadFromDisk(queryPath, filter);
    auto result = classify::retrieveColours(classify::classify(filter, db, (filter.last.size() + 10 * threads - 1) / (10 * threads)), db);
    // write to output
    stream::StreamManager mgr;
    std::ostream &os = mgr[outputFile];
    result.flush(os);
  }
  catch (const exception &e)
  {
    LOG(ERROR) << "Classification failed: " << e.what();
    return 1;
  }
  return 0;
}

// ============================
// DistanceCommand
// ============================
int DistanceCommand::execute()
{
  if (!FilesystemUtils::checkPathExists(databasePath, "Database"))
    return -1;
  graph::StaticGraphDiskFileConfig config(databasePath);
  LOG(INFO) << "Computing distance matrix";
  try
  {
    graph::TraversableGraph db;
    graph::loadSmallBuffers(config.meta, db);
    dist::pairwise::Matrix m = dist::pairwise::countPairsDense(db);
    dist::pairwise::writeToDisk(m, outputFile);
  }
  catch (const exception &e)
  {
    LOG(ERROR) << "Distance computation failed: " << e.what();
    return 1;
  }
  return 0;
}

// ============================
// CountCommand
// ============================
int CountCommand::execute()
{
  if (!FilesystemUtils::checkPathExists(databasePath, "Database"))
    return -1;
  graph::StaticGraphDiskFileConfig config(databasePath);
  LOG(INFO) << "Counting unique k-mers";
  try
  {
    graph::TraversableGraph db;
    graph::loadSmallBuffers(config.meta, db);
    auto data = useWhole
                    ? dist::unique::countUnique(db, dist::bitmask(db.getBaseColourBitWidth()), db.getMaxColour(), 1)
                    : dist::unique::countUnique(db, identity{}, db.getMaxColour(), 1);
    dist::unique::writeToDisk(data, outputFile);
  }
  catch (const exception &e)
  {
    LOG(ERROR) << "Counting failed: " << e.what();
    return 1;
  }
  return 0;
}

// ============================
// CoverageCommand
// ============================
int CoverageCommand::execute()
{
  if (!FilesystemUtils::checkPathExists(inputFile, "Input file"))
    return -1;
  if (!FilesystemUtils::checkPathExists(countFile, "Count file"))
    return -1;
  try
  {
    LOG(INFO) << "parsing unique counts from " << countFile;
    auto uniq = ParseUtils::parseCounts(countFile);
    LOG(INFO) << "accumulating counts from " << inputFile;
    std::unordered_map<size_t, size_t> obs;
    std::vector<size_t> current;
    ParseUtils::parse(inputFile, (size_t)(skipHeader ? 1u : 0u), [&](const std::string &line) -> void
                      {
            ParseUtils::parseNumbersToList(current, line.substr(0, line.find('\t')));
            size_t val = ParseUtils::tabSelectAsSize_t(line, valueColumn);
            for (const auto &k : current) {
                obs[k] += val;
            }
            current.clear(); });
    LOG(INFO) << "writing output to " << outputFile;
    stream::StreamManager mgr;
    ostream &os = mgr[outputFile];
    for (const auto &[k, v] : obs)
    {
      os << k << '\t' << v << '\t' << uniq[k] << '\n';
    }
  }
  catch (const exception &e)
  {
    LOG(ERROR) << "Counting failed: " << e.what();
    return 1;
  }
  return 0;
}

// ============================
// SubcountCommand
// ============================
int SubcountCommand::execute()
{
  if (!FilesystemUtils::checkPathExists(databasePath, "Database"))
    return -1;
  graph::StaticGraphDiskFileConfig config(databasePath);
  LOG(INFO) << "Counting subsamples";
  try
  {
    graph::TraversableGraph db;
    graph::loadSmallBuffers(config.meta, db);
    auto sizes = extractSubsamples(db.getNumGenomes());
    if (sizes.empty())
    {
      LOG(ERROR) << "No valid subsample sizes provided. Please check --subsamples argument.";
      return -1;
    }
    vector<unordered_set<size_t>> randomIds = generateRandomGenomeIdSets(sizes, reps, db.getNumGenomes(), 1);
    unique_ptr<atomic_uint64_t[]> counts = dist::unique::countKmerSubsamples(db, randomIds, 1);
    stream::StreamManager mgr;
    ostream &os = mgr[outputFile];
    for (size_t i = 0; i < randomIds.size(); ++i)
      os << randomIds[i].size() << '\t' << counts[i].load() << '\n';
  }
  catch (const exception &e)
  {
    LOG(ERROR) << "Subcount failed: " << e.what();
    return 1;
  }
  return 0;
}

// ============================
// CompressCommand
// ============================
int CompressCommand::execute()
{
  if (!FilesystemUtils::checkPathExists(databasePath, "Database"))
    return -1;
  graph::StaticGraphDiskFileConfig config(databasePath);
  LOG(INFO) << "Counting co-occurring k-mers";
  try
  {
    graph::TraversableGraph db;
    graph::loadSmallBuffers(config.meta, db);
    auto mat = dist::cooc::countPatterns(db, 1);
    stream::StreamManager mgr;
    std::ostream &os = mgr[outputFile];
    mat.flush(os);
  }
  catch (const exception &e)
  {
    LOG(ERROR) << "Compression failed: " << e.what();
    return 1;
  }
  return 0;
}

// ============================
// JaccardCommand
// ============================
int JaccardCommand::execute()
{
  if (!FilesystemUtils::checkPathExists(databasePath, "Database"))
    return -1;
  if (!FilesystemUtils::checkPathExists(nwkFile, "Phylogeny tree"))
    return -1;
  if (!FilesystemUtils::checkPathExists(countFile, "Count file"))
    return -1;
  if (!FilesystemUtils::checkPathExists(distFile, "Distance file"))
    return -1;
  graph::StaticGraphDiskFileConfig config(databasePath);
  LOG(INFO) << "Computing Jaccard statistics";
  try
  {
    graph::TraversableGraph db;
    graph::loadSmallBuffers(config.meta, db);
    LOG(INFO) << "loading genome manifest";
    Metadata::GenomeManifest genomeIds = Metadata::loadGenomeManifest(config.genomesids);
    LOG(INFO) << "parsing counts";
    auto uniqKmers = ParseUtils::parseCounts(countFile);
    LOG(INFO) << "loading phylogeny";
    PhylogenyUtils::Phylogeny phy = PhylogenyUtils::loadAndReduce(nwkFile, genomeIds);
    double z = stats::detail::normalPPF(1.0 - alpha / 2.0);
    stream::StreamManager mgr;
    LOG(INFO) << "making Jaccard sampler";
    auto os = stats::makeJaccardSampler(take, z, mgr[outputFile]);
    LOG(INFO) << "writing output";
    ParseUtils::parse(distFile, [&](const std::string &line)
                      {
            stats::JaccardResult rec;
            rec.gid1 = ParseUtils::tabSelectAsSize_t(line, 0);
            rec.gid2 = ParseUtils::tabSelectAsSize_t(line, 1);
            rec.k = db.k + 1;
            rec.shared = ParseUtils::tabSelectAsSize_t(line, 2);
            const GenomeToken &info1 = genomeIds.ids.at(rec.gid1), &info2 = genomeIds.ids.at(rec.gid2);
            rec.c1 = uniqKmers.at(rec.gid1);
            rec.c2 = uniqKmers.at(rec.gid2);
            rec.l1 = info1.sequenceLength;
            rec.l2 = info2.sequenceLength;
            double shd = static_cast<double>(rec.shared), u1 = static_cast<double>(rec.c1), u2 = static_cast<double>(rec.c2);
            rec.jcd = shd / (u1 + u2 - shd);
            rec.lcaName = PhylogenyUtils::getNearestNamedLca(info1.name, info2.name, phy);
            os->write(rec); });
  }
  catch (const exception &e)
  {
    LOG(ERROR) << "Jaccard computation failed: " << e.what();
    return 1;
  }
  return 0;
}

// ============================
// LcaCommand
// ============================
int LcaCommand::execute()
{
  if (!FilesystemUtils::checkPathExists(databasePath, "Database"))
    return -1;
  if (!FilesystemUtils::checkPathExists(clusterFile, "Cluster file"))
    return -1;
  if (!FilesystemUtils::checkPathExists(countFile, "Count file"))
    return -1;
  if (!FilesystemUtils::checkPathExists(nwkFile, "Phylogeny tree"))
    return -1;
  graph::StaticGraphDiskFileConfig config(databasePath);
  try
  {
    LOG(INFO) << "Loading graph metadata";
    // load metadata
    graph::TraversableGraph g;
    graph::loadSmallBuffers(config.meta, g);
    uint8_t baseIdWidth = g.getBaseColourBitWidth();
    dist::bitmask msk(baseIdWidth);
    LOG(INFO) << "Loading Genome IDs";
    Metadata::GenomeManifest genomeIds = Metadata::loadGenomeManifest(config.genomesids);
    LOG(INFO) << "Loading Annotation IDs";
    Metadata::AnnotationManifest annotIds = Metadata::loadAnnotationManifest(config.annotids);
    PhylogenyUtils::Phylogeny phy = PhylogenyUtils::loadAndReduce(nwkFile, genomeIds);
    ClusterUtils::Clustering clusters(baseIdWidth);
    LOG(INFO) << "Collecting Annotation ID seeds";
    clusters.loadAllSeeds(genomeIds, annotIds);
    LOG(INFO) << "Loading seed cluster assignments";
    clusters.loadClusters(clusterFile);
    // process counts file
    LOG(INFO) << "Computing LCA conservation";
    stream::StreamManager mgr;
    std::ostream &os = mgr[outputFile];
    if (printColours)
      os << "annot_ids\t";
    os << "group_size\tsignature_size\tnum_coding\tnum_families\tfamily_name\tnum_in_tree\ttree_clade_size\tlca\n";
    ParseUtils::Lca1::CompressItem result;
    ParseUtils::parse(countFile, [&](const std::string &line) -> void
                      {
            // collect cluster data
            result.parseColours(line, clusters);
            if (result.numColours == 0) return;
            // check family diversity
            if (result.numCoding > 0) { // all non-coding are allowed to pass by default
                float famDiv = (static_cast<float>(result.numClusters()) - 1.0f)
                    / (static_cast<float>(result.numColours) - 1.0f);
                if ((famDiv < minFamDiv) || (famDiv > maxFamDiv))
                    return;
            }
            // collect conservation data
            result.collectPhylogenyNodes(msk, genomeIds, phy);
            if (result.cladeSize == 0) return;
            // check conservation
            float phyCons = static_cast<float>(result.numGenomes())
                / static_cast<float>(result.cladeSize);
            if ((phyCons < minPhyCons) || (phyCons > maxPhyCons))
                return;
            // print results
            if (printColours) {
                std::string_view sv = line;
                sv = sv.substr(0, sv.find('\t'));
                os << sv << '\t';
            }
            os << result; });
  }
  catch (const exception &e)
  {
    LOG(ERROR) << "LCA computation failed: " << e.what();
    return 1;
  }
  return 0;
}

std::ostream &operator<<(std::ostream &os, const ParseUtils::Lca1::CompressItem &obj)
{
  return obj.print(os);
}

// ============================
// GlcaCommand
// ============================
int GlcaCommand::execute()
{
  if (!FilesystemUtils::checkPathExists(databasePath, "Database"))
    return -1;
  if (!FilesystemUtils::checkPathExists(clusterFile, "Cluster file"))
    return -1;
  if (!FilesystemUtils::checkPathExists(nwkFile, "Phylogeny tree"))
    return -1;
  graph::StaticGraphDiskFileConfig config(databasePath);
  LOG(INFO) << "Computing gene conservation";
  try
  {
    graph::TraversableGraph g;
    graph::loadSmallBuffers(config.meta, g);
    Metadata::AccessionManifest accessions = Metadata::collectAccessions(config.annotids, g.getBaseColourBitWidth());
    Metadata::GenomeManifest genomeIds = Metadata::loadGenomeManifest(config.genomesids);
    PhylogenyUtils::Phylogeny phy = PhylogenyUtils::loadAndReduce(nwkFile, genomeIds);
    ClusterUtils::GenomeClusters clusters = ClusterUtils::groupByGenomeId(clusterFile, accessions, genomeIds, phy);
    stream::StreamManager mgr;
    std::ostream &os = mgr[outputFile];
    for (auto &[seed, nodes] : clusters)
    {
      if (nodes.empty())
        continue;
      ParseUtils::unique(nodes);
      size_t lcaId, numChildren;
      if (nodes.size() == 1)
      {
        lcaId = nodes[0];
        numChildren = 1;
      }
      else
      {
        lcaId = phy.lca(nodes);
        numChildren = phy.tree.countSubleaves(lcaId);
      }
      os << seed << '\t'
         << nodes.size() << '\t'
         << numChildren << '\t'
         << PhylogenyUtils::getNearestNamedAncestor(lcaId, phy) << '\n';
    }
  }
  catch (const exception &e)
  {
    LOG(ERROR) << "GLCA computation failed: " << e.what();
    return 1;
  }
  return 0;
}

// ============================
// SamdistCommand
// ============================
int SamdistCommand::execute()
{
  if (!FilesystemUtils::checkPathExists(databasePath, "Database"))
    return -1;
  if (!FilesystemUtils::checkPathExists(countFile, "Count file"))
    return -1;
  if (!FilesystemUtils::checkPathExists(nwkFile, "Phylogeny tree"))
    return -1;
  graph::StaticGraphDiskFileConfig config(databasePath);
  LOG(INFO) << "Computing gene distances";
  try
  {
    graph::TraversableGraph db;
    graph::loadSmallBuffers(config.meta, db);
    uint8_t baseIdWidth = db.getBaseColourBitWidth();
    dist::pairwise::ColourSampler mat = dist::pairwise::samplePairsDense(db);
    Metadata::GenomeManifest genomeIds = Metadata::loadGenomeManifest(config.genomesids);
    Metadata::AnnotationManifest annotIds = Metadata::loadAnnotationManifest(config.annotids);
    PhylogenyUtils::Phylogeny phy = PhylogenyUtils::loadAndReduce(nwkFile, genomeIds);
    ClusterUtils::Clustering clusters(baseIdWidth);
    clusters.loadSeeds(mat.gatherKeys(baseIdWidth), genomeIds, annotIds);
    auto uniqKmers = ParseUtils::parseCounts(countFile);
    {
      stream::StreamManager mgr;
      std::ostream &os = mgr[outputFile];
      StatCalculator::SummariseColourPairData Op(baseIdWidth, db.k + 1, uniqKmers, genomeIds, phy, clusters, os);
      oneapi::tbb::enumerable_thread_specific tlocal(Op);
      mat.pforEach(StatCalculator::SummariseColourPairData::capacity, tlocal);
    }
  }
  catch (const exception &e)
  {
    LOG(ERROR) << "Samdist computation failed: " << e.what();
    return 1;
  }
  return 0;
}

// ============================
// QueryCommand
// ============================
int QueryCommand::execute()
{
  if (!FilesystemUtils::checkPathExists(queryFile, "Query"))
    return -1;
  if (!FilesystemUtils::checkPathExists(databasePath, "Database"))
    return -1;
  if (!FilesystemUtils::checkPathExists(countFile, "Count file"))
    return -1;
  if (!FilesystemUtils::checkPathExists(nwkFile, "Phylogeny tree"))
    return -1;
  graph::StaticGraphDiskFileConfig config(databasePath);
  try
  {
    // load small metadata
    graph::TraversableGraph db;
    graph::loadSmallBuffers(config.meta, db);
    Metadata::GenomeManifest genomeIds = Metadata::loadGenomeManifest(config.genomesids);
    PhylogenyUtils::Phylogeny phy = PhylogenyUtils::loadAndReduce(nwkFile, genomeIds);
    // perform search
    dist::query::QuerySet queries = QueryUtils::parseQueries(queryFile, genomeIds);
    dist::query::ConcurrentBuffers results = dist::query::searchGenomes(queries, db);
    // load larger metadata
    Metadata::AnnotationManifest annotIds = Metadata::loadAnnotationManifest(config.annotids);
    auto uniqKmers = ParseUtils::parseCounts(countFile);
    // load cluster data for search results
    ClusterUtils::Clustering clusters(db.getBaseColourBitWidth());
    clusters.loadSeeds(results.gatherKeys(), genomeIds, annotIds);
    // write results to output
    {
      dist::bitmask msk(db.getBaseColourBitWidth());
      stream::StreamManager mgr;
      std::ostream &os = mgr[outputFile];
      oneapi::tbb::parallel_for(
          (size_t)0, results.size(), (size_t)1,
          [&](size_t bufferIndex) -> void
          {
            StatCalculator::SummariseColourPairData Summarise(db.getBaseColourBitWidth(), db.k + 1, uniqKmers, genomeIds, phy, clusters, os);
            for (const auto &[kp, v] : results[bufferIndex].data)
              Summarise(msk(kp.first), msk(kp.second), kp.first, kp.second, v.load());
          });
    }
  }
  catch (const exception &e)
  {
    LOG(ERROR) << "Query computation failed: " << e.what();
    return 1;
  }
  return 0;
}

// ============================
// KmersCommand
// ============================
int KmersCommand::execute()
{
  if (!FilesystemUtils::checkPathExists(coloursFile, "Colours"))
    return -1;
  if (!FilesystemUtils::checkPathExists(databasePath, "Database"))
    return -1;
  graph::StaticGraphDiskFileConfig config(databasePath);
  try
  {
    // load graph into RAM
    LOG(INFO) << "loading graph into RAM";
    graph::TraversableGraph db;
    graph::loadFromDisk(config, db);
    // load colour profiles from disk
    LOG(INFO) << "reading colour profiles from first column of " << coloursFile;
    auto colourProfiles = ParseUtils::parseColourProfiles(coloursFile, db.getColourBitWidth(), skipHeader);
    // retrieve colour profiles
    LOG(INFO) << "starting graph buffer traversal";
    stream::StreamManager mgr;
    std::ostream &os = mgr[outputFile];
    dist::cooc::printColours(db, colourProfiles, os);
  }
  catch (const exception &e)
  {
    LOG(ERROR) << "Kmers computation failed: " << e.what();
    return 1;
  }
  return 0;
}

// ============================
// Common workflow implementations
// ============================
namespace LoggingUtils
{
  void setupLogging(const string &logFile, const char *execName)
  {
    if (logFile.empty())
    {
      FLAGS_logtostderr = 1;
    }
    else
    {
      for (auto severity : {google::INFO, google::WARNING, google::ERROR, google::FATAL})
      {
        google::SetLogDestination(severity, logFile.c_str());
        google::SetLogSymlink(severity, "");
      }
      FLAGS_alsologtostderr = true;
    }
    google::InitGoogleLogging(execName);
  }
}

namespace ThreadUtils
{
  void configureThreads(size_t threads)
  {
    static oneapi::tbb::global_control global_limit(oneapi::tbb::global_control::max_allowed_parallelism, threads);
  }
}

namespace FilesystemUtils
{
  bool checkPathExists(const string &path, const string &description)
  {
    if (!fs::exists(path))
    {
      LOG(ERROR) << description << " not found at: " << path << ". Please check the path and permissions.";
      return false;
    }
    return true;
  }

  bool ensureDirectoryExists(const string &path, bool allowOverwrite)
  {
    fs::path fpath = path;
    if (!fs::exists(fpath))
      fs::create_directory(fpath);
    else if (!allowOverwrite)
    {
      LOG(ERROR) << "directory already exists: " << path;
      return false;
    }
    return true;
  }
}

namespace ValidationUtils
{
  bool validateIndexingFlags(bool useBulk, bool suffixSize, long double resourceLimit)
  {
    if (useBulk && suffixSize > 0)
    {
      LOG(ERROR) << "flags --bulk and --suffix-size are incompatible, please only supply one";
      return false;
    }
    else if (useBulk && resourceLimit > 0.0l)
    {
      LOG(ERROR) << "flags --bulk and --ram-gb are incompatible, please only supply one";
      return false;
    }
    else if (suffixSize > 0 && resourceLimit > 0.0l)
    {
      LOG(ERROR) << "flags --suffix-size and --ram-gb are incompatible, please only supply one";
      return false;
    }
    return true;
  }
}

namespace ParseUtils
{
  size_t tab(const std::string &s, size_t t)
  {
    size_t pos = 0;
    for (size_t occs = 0; pos < s.size() && occs < t; ++pos)
      if (s[pos] == '\t')
        ++occs;
    if (pos == s.size())
      throw std::runtime_error(std::string{"Index out of range: "} + std::to_string(t));
    return pos;
  }

  size_t tabFind(std::string const &s, const char *qry)
  {
    size_t index = s.find(qry);
    if (index == std::string::npos)
      throw std::runtime_error(std::string{"Could not find header column"} + qry);
    return std::count(s.cbegin(), s.cbegin() + index, '\t');
  }

  std::pair<size_t, size_t> tabSelect(std::string const &s, size_t t)
  {
    size_t begin = tab(s, t), end = s.find('\t', begin + 1);
    if (end == std::string::npos)
      end = s.size();
    return std::make_pair(begin, end);
  }

  size_t tabSelectAsSize_t(std::string const &s, size_t t)
  {
    auto rng = tabSelect(s, t);
    size_t val;
    try
    {
      val = std::stoull(s.substr(rng.first, rng.second - rng.first));
    }
    catch (const std::invalid_argument &e)
    {
      throw std::runtime_error("Could not size_t from " + s.substr(rng.first, rng.second - rng.first));
    }
    return val;
  }

  std::string tabExtractString(const std::string &s, size_t t)
  {
    auto rng = tabSelect(s, t);
    return s.substr(rng.first, rng.second - rng.first);
  }

  std::unordered_map<size_t, size_t> parseCounts(const std::string &f)
  {
    std::unordered_map<size_t, size_t> data;
    parse(f, [&](const std::string &line)
          { data[tabSelectAsSize_t(line, 0)] = tabSelectAsSize_t(line, 1); });
    return data;
  }

  size_t Lca1::CompressItem::numClusters() const noexcept
  {
    return clusters.size();
  }

  size_t Lca1::CompressItem::numGenomes() const noexcept
  {
    return values.size();
  }

  void Lca1::CompressItem::parseColours(const std::string &line, const ClusterUtils::Clustering &lookup)
  {
    itemValue = tabSelectAsSize_t(line, 1);
    values.clear();
    parseNumbersToList(values, line.substr(0, line.find('\t')));
    clusters.clear();
    numCoding = 0;
    numColours = values.size();
    for (const uint64_t &aid : values)
    {
      if (const std::string *rep = lookup.rep(aid); rep != nullptr)
      {
        ++numCoding;
        for (const std::string *&otherRep : clusters)
          if (rep == otherRep)
            goto finishLookup;
        clusters.push_back(rep);
      finishLookup:
      }
    }
  }

  void Lca1::CompressItem::collectPhylogenyNodes(const dist::bitmask &msk, const Metadata::GenomeManifest &genomeIds, const PhylogenyUtils::Phylogeny &phy)
  {
    PhylogenyUtils::colourToNodeId(values, msk, genomeIds, phy);
    if (values.empty())
    {
      cladeSize = 0; // used for filtering
      lcaName.clear();
    }
    else
    {
      unique(values);
      size_t lcaId = phy.lca(values);
      cladeSize = phy.tree.countSubleaves(lcaId);
      lcaName = PhylogenyUtils::getNearestNamedAncestor(lcaId, phy);
    }
  }

  std::ostream &Lca1::CompressItem::print(std::ostream &os) const
  {
    os << numColours << '\t'
       << itemValue << '\t'
       << numCoding << '\t'
       << numClusters() << '\t';
    if (numClusters() == 1)
      os << *clusters[0];
    os << '\t';
    if (values.empty())
      os << "0\t0\t\n";
    else
    {
      os << values.size() << '\t'
         << cladeSize << '\t'
         << lcaName << '\n';
    }
    return os;
  }

  dist::cooc::RankMap parseColourProfiles(const std::string &colourFile, size_t bitWidth, bool skipHeader)
  {
    dist::cooc::RankMap map(bitWidth);
    dist::ColourVector tmp;
    size_t i = 0;
    parse(colourFile, (size_t)(skipHeader ? 1u : 0u), [&](const std::string &line) -> void
          {
            parseNumbersToList(tmp, line.substr(0, line.find_first_of("\t\n ")));
            map.set(tmp, i++);
            tmp.clear(); });
    return map;
  }

}

namespace GraphBuilder
{
  IndexParams::IndexParams(const std::string &queryFile_, const std::string &outputFolder_, uint8_t kmerSize_, size_t threads_, long double resourceLimit_, bool useBulk_, uint8_t suffixSize_)
      : queryFile(queryFile_), outputFolder(outputFolder_), kmerSize(kmerSize_), threads(threads_), resourceLimit(resourceLimit_), useBulk(useBulk_), suffixSize(suffixSize_), edgeBound(), parsedBytes(), graphBytes(), roundedResourceLimit()
  {
  }

  bool IndexParams::validateBytesLimit()
  {
    if (resourceLimit)
    {
      LOG(INFO) << "Resources used: sequence data (" << static_cast<long double>(parsedBytes) / 1.0e9l << " Gb), predicted graph size (" << static_cast<long double>(graphBytes) / 1.0e9l << " Gb)";
      size_t totalBytesUsed = parsedBytes + graphBytes;
      roundedResourceLimit = static_cast<size_t>(resourceLimit * 1.0E9l);
      if (roundedResourceLimit < totalBytesUsed)
      {
        LOG(ERROR) << "resource limit is less than est. graph size: " << static_cast<long double>(totalBytesUsed) / 1.0e9l << "Gb";
        return false;
      }
      roundedResourceLimit -= totalBytesUsed;
    }
    return true;
  }

  bool buildIndex(IndexParams &params)
  {
    auto genomeFiles = readFilePaths(params.queryFile.data());
    graph::GraphDiskFileConfig config(params.outputFolder);
    FeatureIdGenerator<GenomeStatsSummary, Gff3Record> encoder(genomeFiles.size());
    Dna4GenomeVector genomes = loadGenomes(genomeFiles, encoder, params.kmerSize + 1u);
    encoder.writeFeatureNamesTo(config.genomesids.string());
    encoder.clearFeatureNames();
    encoder.writeFeaturesTo(config.annotids.string());
    encoder.clearFeatures();
    try
    {
      graph::WriteableGraph garr = buildIndexInMemory(params, genomes, config);
      genomes.reset_memory();
      graph::makeStaticGraphOnDisk(garr, params.outputFolder);
    }
    catch (const ArgumentInvalidationException &)
    {
      return false;
    }
    return true;
  }

  bool buildFilter(IndexParams &params)
  {
    fs::path db_dir = params.outputFolder;
    if (!fs::exists(db_dir))
      fs::create_directory(db_dir);
    string data(istreambuf_iterator<char>(cin), {});
    istringstream din(data);
    Dna4Genome genome{};
    parseFastaStream(genome, din, 0, params.kmerSize + 1u);
    if (!genome.length())
    {
      LOG(ERROR) << "empty input FASTA.";
      return false;
    }
    ChunkedDna4Genome xgenome(move(genome), params.threads, params.kmerSize);
    try
    {
      graph::BaseGraph garr = buildFilterInMemory(params, xgenome);
      xgenome.reset_memory();
      graph::makeStaticGraphOnDisk(garr, params.outputFolder);
    }
    catch (const ArgumentInvalidationException &)
    {
      return false;
    }
    return true;
  }
}

namespace ClusterUtils
{
  Clustering::Clustering(uint8_t idWidth_) : baseIdWidth(idWidth_) {}

  size_t Clustering::StrHash::operator()(const char *s) const
  {
    return hash_type{}(s);
  }

  size_t Clustering::StrHash::operator()(std::string_view s) const
  {
    return hash_type{}(s);
  }

  size_t Clustering::StrHash::operator()(const std::string &s) const
  {
    return hash_type{}(s);
  }

  void Clustering::insert(uint64_t qry, std::string &&seed)
  {
    auto [it, _] = seeds.emplace(std::move(seed));
    lookup[qry] = std::to_address(it);
  }

  void Clustering::assignIfExist(std::string_view rep, std::string_view seed)
  {
    auto it = seeds.find(seed);
    if (it != seeds.end())
    {
      const std::string *from = std::to_address(it);
      auto [repit, _] = seeds.emplace(rep);
      reps[from] = std::to_address(repit);
    }
  }

  const std::string *Clustering::seed(uint64_t qry) const
  {
    if (auto it = lookup.find(qry); it != lookup.end())
      return it->second;
    else
      return nullptr;
  }

  const std::string *Clustering::rep(uint64_t qry) const
  {
    const std::string *s = seed(qry);
    if (s)
    {
      if (auto it = reps.find(s); it != reps.end())
        return it->second;
      else
        return nullptr;
    }
    else
      return nullptr;
  }

  void Clustering::insertGenomeSeeds(const Metadata::GenomeManifest &genomeIds)
  {
    for (const auto &[id, tok] : genomeIds.ids)
    {
      insert(id, std::string{tok.name});
      assignIfExist(tok.name, tok.name);
    }
  }

  void Clustering::loadAllSeeds(const Metadata::GenomeManifest &genomeIds, const Metadata::AnnotationManifest &annotIds)
  {
    insertGenomeSeeds(genomeIds);
    auto requests = Metadata::retrieveAllSources(baseIdWidth, genomeIds, annotIds);
    retrieveRequests(requests);
  }

  void Clustering::loadSeeds(const std::unordered_set<size_t> &ids, const Metadata::GenomeManifest &genomeIds, const Metadata::AnnotationManifest &annotIds)
  {
    insertGenomeSeeds(genomeIds);
    auto requests = Metadata::retrieveSources(ids, baseIdWidth, genomeIds, annotIds);
    retrieveRequests(requests);
  }

  void Clustering::retrieveRequests(const std::unordered_map<std::string, Metadata::GffRequestSet> &requests)
  {
    std::mutex mtx_;
    oneapi::tbb::parallel_for(
        (size_t)0, requests.bucket_count(), [&](size_t bucketIndex) -> void
        {
                std::vector<std::pair<uint64_t,std::string>> threadCache; // result cache
                for (auto r = requests.begin(bucketIndex); r != requests.end(bucketIndex); ++r) {
                    const Metadata::GffRequestSet &thisRequests = r->second;
                    // process file
                    zstr::ifstream gffStream(r->first);
                    std::string line;
                    while (std::getline(gffStream, line)) {
                        if (line.starts_with('#')) {
                            if (line == "##FASTA") break;
                        } else {
                            Gff3Record rec(line);
                            if (auto thisRequest = thisRequests.find(rec); thisRequest != thisRequests.cend()) {
                                // fetch seed name
                                std::string seed = rec.ref;
                                if (size_t parentPos = line.find("Parent=", line.find_last_of('\t')+1); parentPos != std::string::npos) {
                                    seed += '-';
                                    size_t endPos = line.find(';', parentPos + 7);
                                    if (endPos == std::string::npos) endPos = line.size();
                                    seed.append(line, parentPos + 7, endPos - (parentPos + 7));
                                }
                                // add result to cache
                                threadCache.emplace_back(thisRequest->second, std::move(seed));
                            }
                        }
                    }
                    // flush cache under lock protection
                    {
                        std::lock_guard lock(mtx_);
                        for (auto &[qry, seed] : threadCache) {
                            insert(qry, std::move(seed));
                        }
                    }
                    threadCache.clear();
                } });
  }

  void ClusterUtils::Clustering::loadClusters(const std::string &file)
  {
    ParseUtils::parse(file, [&](std::string_view line)
                      {
            size_t sep = line.find('\t');
            if (sep != std::string::npos) assignIfExist(line.substr(0, sep), line.substr(sep+1, line.size()-(sep+1)));
            else throw std::runtime_error("Unformatted line: " + std::string(line)); });
  }

  std::string commonCluster(const Clustering &clusters, uint64_t id1, uint64_t id2)
  {
    const std::string *rep1 = clusters.rep(id2), *rep2 = clusters.rep(id2);
    if (!rep1)
      throw std::runtime_error("Could not find cluster for id " + std::to_string(id1));
    else if (!rep2)
      throw std::runtime_error("Could not find cluster for id " + std::to_string(id2));
    else if (*rep1 != *rep2)
      return "NA";
    else
      return *rep1;
  }

  GenomeClusters groupByGenomeId(const std::string &clusterFile, const Metadata::AccessionManifest &acc, const Metadata::GenomeManifest &genomeIds, const PhylogenyUtils::Phylogeny &phy)
  {
    GenomeClusters clusters;
    ParseUtils::parse(clusterFile, [&](const std::string &line) -> void
                      {
            size_t sep = line.find('\t');
            if (sep == std::string::npos) return;
            auto &cluster = clusters[line.substr(0, sep)];
            size_t end = line.find("-", sep + 1);
            if (end == std::string::npos)
                throw std::runtime_error("Expected cluster name format ACCESSION-...-GENE, could not parse: " + line);
            auto gidp = acc.find(line.substr(sep+1, end-(sep+1)));
            if (gidp == acc.end())
                throw std::runtime_error("Unrecognised accession " + line.substr(sep+1, end-(sep+1)));
            auto gnamep = genomeIds.ids.find(gidp->second);
            if (gnamep == genomeIds.ids.end())
                throw std::runtime_error("Unrecognised genome ID " + std::to_string(gidp->second));
            if (auto node = phy.tree.id(gnamep->second.name); node != phy.tree.end())
                cluster.push_back(node->second); });
    return clusters;
  }
}

namespace Metadata
{
  GenomeManifest loadGenomeManifest(const std::string &file)
  {
    GenomeManifest manifest;
    bool inHeader = true;
    size_t idcol = 0, namecol = 1, srccol = 2, lencol = 6;
    ParseUtils::parse(file, [&](const std::string &line)
                      {
            if (inHeader) {
                idcol = ParseUtils::tabFind(line, "sequenceID");
                namecol = ParseUtils::tabFind(line, "genomeName");
                srccol = ParseUtils::tabFind(line, "sourceFile");
                lencol = ParseUtils::tabFind(line, "sequenceLength");
                inHeader = false;
            } else {
                manifest.ids.try_emplace(
                    ParseUtils::tabSelectAsSize_t(line, idcol),
                    ParseUtils::tabExtractString(line, namecol),
                    ParseUtils::tabSelectAsSize_t(line, lencol),
                    ParseUtils::tabExtractString(line, srccol)
                );
            } });
    return manifest;
  }

  AnnotationManifest loadAnnotationManifest(const std::string &file)
  {
    AnnotationManifest manifest;
    bool inHeader = true;
    size_t idcol = 0, accncol = 1, bgcol = 2, edcol = 3, strcol = 4;
    ParseUtils::parse(file, [&](const std::string &line)
                      {
            if (inHeader) {
                idcol = ParseUtils::tabFind(line, "annotID");
                accncol = ParseUtils::tabFind(line, "accession");
                bgcol = ParseUtils::tabFind(line, "begin");
                edcol = ParseUtils::tabFind(line, "end");
                strcol = ParseUtils::tabFind(line, "strand");
                inHeader = false;
            } else {
                manifest.ids.try_emplace(
                    ParseUtils::tabSelectAsSize_t(line, idcol),
                    ParseUtils::tabExtractString(line, accncol),
                    ParseUtils::tabSelectAsSize_t(line, bgcol),
                    ParseUtils::tabSelectAsSize_t(line, edcol),
                    line[ParseUtils::tab(line, strcol)]
                );
            } });
    return manifest;
  }

  AccessionManifest collectAccessions(const std::string &file, uint8_t idWidth)
  {
    AccessionManifest accessions;
    dist::bitmask msk(idWidth);
    ParseUtils::parse(file, [&](const std::string &line)
                      {
            size_t annotId = ParseUtils::tabSelectAsSize_t(line, 0);
            auto [accnBegin, accnEnd] = ParseUtils::tabSelect(line, 1);
            accessions.try_emplace(line.substr(accnBegin, accnEnd - accnBegin), msk(annotId)); });
    return accessions;
  }

  std::unordered_map<std::string, GffRequestSet> retrieveAllSources(uint8_t baseIdWidth, const GenomeManifest &genomeIds, const AnnotationManifest &annotIds)
  {
    std::unordered_map<std::string, GffRequestSet> requests;
    dist::bitmask msk(baseIdWidth);
    for (const auto &[id, annot] : annotIds.ids)
    {
      if (id == msk(id))
        continue;
      auto grec = genomeIds.ids.find(msk(id));
      if (grec == genomeIds.ids.cend())
        throw std::runtime_error("Unrecognised genome ID " + std::to_string(msk(id)));
      requests[grec->second.source][annot] = id;
    }
    return requests;
  }
}

namespace PhylogenyUtils
{
  Phylogeny::Phylogeny(const std::string &nwk_) : tree(nwk_.c_str()), lcaIndex(tree) {}
  Phylogeny::Phylogeny(maki_tree::Tree &&tree_) : tree(std::move(tree_)), lcaIndex(tree) {}

  Phylogeny loadAndReduce(const std::string &nwkFile, const Metadata::GenomeManifest &genomes)
  {
    maki_tree::Tree temp;
    temp.load(nwkFile.c_str());
    size_t missing = 0;
    std::unordered_set<std::string> currentGenomes;
    for (const auto &[_, tk] : genomes.ids)
    {
      if (temp.id(tk.name) != temp.end())
        currentGenomes.insert(tk.name);
      else
        ++missing;
    }
    if (missing)
      LOG(ERROR) << missing << " genomes could not be found in provided phylogeny";
    return Phylogeny(temp.reduce(currentGenomes));
  }

  std::string getNearestNamedLca(const std::string &n1, const std::string &n2, const Phylogeny &phy)
  {
    auto it1 = phy.tree.id(n1), it2 = phy.tree.id(n2);
    if (it1->first != n1 || it2->first != n2)
      return std::string{"NA"};
    size_t lcaId = phy.lca(it1->second, it2->second);
    return getNearestNamedAncestor(lcaId, phy);
  }

  std::string getNearestNamedAncestor(size_t nodeId, const Phylogeny &phy)
  {
    auto resultIt = phy.tree.name(nodeId);
    size_t nameAnchor = resultIt->first.rfind("__");
    while (nodeId != phy.tree.root() && nameAnchor == std::string::npos)
    {
      nodeId = phy.tree.parent(nodeId);
      resultIt = phy.tree.name(nodeId);
      nameAnchor = resultIt->first.rfind("__");
    }
    if (nameAnchor == std::string::npos)
      return resultIt->first; // root
    else
    {
      const std::string &lcaName = resultIt->first;
      size_t start = nameAnchor - 1, end = lcaName.find(';', nameAnchor + 2);
      if (end == std::string::npos)
        end = lcaName.size();
      return lcaName.substr(start, end - start);
    }
  }

  void colourToNodeId(std::vector<size_t> &data, const dist::bitmask &msk, const Metadata::GenomeManifest &genomeIds, const Phylogeny &phy)
  {
    auto out = data.begin();
    for (size_t x : data)
    {
      try
      {
        if (auto it = phy.tree.id(genomeIds.ids.at(msk(x)).name); it != phy.tree.end())
          *out++ = it->second;
      }
      catch (const std::out_of_range &e)
      {
        throw std::runtime_error("Unrecognised colour " + std::to_string(msk(x)));
      }
    }
    data.erase(out, data.end());
  }
}

namespace StatCalculator
{
  SummariseColourPairData::SummariseColourPairData(size_t baseIdWidth_, size_t k_, const std::unordered_map<size_t, size_t> &uniqKmers_, const Metadata::GenomeManifest &genomeIds_, const PhylogenyUtils::Phylogeny &phy_, const ClusterUtils::Clustering &clusters_, std::ostream &os_)
      : records(), size(0), baseIdWidth(baseIdWidth_), k(k_), uniqKmers(uniqKmers_), genomeIds(genomeIds_), phy(phy_), clusters(clusters_), os(os_)
  {
  }

  SummariseColourPairData::SummariseColourPairData(const StatCalculator::SummariseColourPairData &rhs)
      : records(), size(0), baseIdWidth(rhs.baseIdWidth), k(rhs.k), uniqKmers(rhs.uniqKmers), genomeIds(rhs.genomeIds), phy(rhs.phy), clusters(rhs.clusters), os(rhs.os)
  {
  }

  SummariseColourPairData::~SummariseColourPairData()
  {
    flush();
  }

  void SummariseColourPairData::operator()(size_t gid1, size_t gid2, size_t fid1, size_t fid2, size_t shared)
  {
    if (full())
      flush();
    assert(size < capacity);
    stats::SamdistResult &rec = records[size++];
    rec.fid1 = fid1;
    rec.fid2 = fid2;
    rec.k = k;
    rec.c1 = rec.l1 = uniqKmers.at(fid1);
    rec.c2 = rec.l2 = uniqKmers.at(fid2);
    if (shared > 0)
    {
      rec.shared = shared;
      double shd = static_cast<double>(rec.shared), u1 = static_cast<double>(rec.c1), u2 = static_cast<double>(rec.c2);
      rec.jcd = shd / (u1 + u2 - shd);
    }
    rec.lcaName = PhylogenyUtils::getNearestNamedLca(genomeIds.ids.at(gid1).name, genomeIds.ids.at(gid2).name, phy);
    rec.seed1 = clusters.seed(fid1);
    rec.seed2 = clusters.seed(fid2);
  }

  void SummariseColourPairData::operator()(size_t gid1, size_t gid2, const dist::pairwise::ColourSampler::KeyValuePair &kv)
  {
    size_t fid1 = gid1 | ((kv.key & 0xFFFFFFFFULL) << baseIdWidth), fid2 = gid2 | ((kv.key & ~0xFFFFFFFFULL) >> (32u - baseIdWidth));
    return operator()(gid1, gid2, fid1, fid2, kv.value);
  }

  bool SummariseColourPairData::full() const
  {
    return size >= capacity;
  }

  void SummariseColourPairData::flush()
  {
    if (size)
    {
      std::osyncstream dest(os);
      for (size_t i = 0; i < size; ++i)
      {
        dest << records[i] << '\n';
        records[i].clear();
      }
    }
    size = 0;
  }
}

namespace ReadsUtils
{
  vector<reads::DataFilePair> pairReads(const string &readsFile)
  {
    auto inputFiles = readFilePaths(readsFile.data());
    return reads::pairInputFiles(inputFiles);
  }

  vector<reads::ReadChunks> chunkReads(const reads::DataFilePair &fp, size_t minReadLength, size_t minFragLength, size_t threads)
  {
    return reads::chunk(reads::parse::parsePairedFastq(fp, minReadLength, minFragLength, threads), threads);
  }
}

namespace QueryUtils
{
  dist::query::QuerySet parseQueries(const std::string &queryFile, const Metadata::GenomeManifest &genomeIds)
  {
    dist::query::QuerySet queries;
    std::unordered_set<std::string> requestedNames;
    // collect names from file
    ParseUtils::parse(queryFile, [&](const std::string &line)
                      { requestedNames.emplace(ParseUtils::tabExtractString(line, 0)); });
    // collect ids
    for (const auto &[gid, rec] : genomeIds.ids)
    {
      if (auto it = requestedNames.find(rec.name); it != requestedNames.end())
      {
        queries.insert(gid);
        requestedNames.erase(it);
      }
    }
    // check for missing ids
    if (requestedNames.size() > 0)
    {
      LOG(ERROR) << "could not identify " << requestedNames.size() << " requested genomes";
    }
    return queries;
  }
}
