// Offline all-node sequence + distinct logical GBWT path count index.
// Standard C++17; every standard header is included explicitly (newer GCC / libc++ no longer
// include <cstdint> or <string> transitively). Compiled by graph_index.py `compile`.
#include <gbwtgraph/gbz.h>
#include <nlohmann/json.hpp>
#include <sqlite3.h>
#include <chrono>
#include <cstdint>
#include <exception>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>
#include <stdexcept>
using Clock=std::chrono::steady_clock;
static double elapsed(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}
static void sql(sqlite3* db,const char* s){char* e=nullptr;if(sqlite3_exec(db,s,nullptr,nullptr,&e)!=SQLITE_OK){std::string m=e?e:"SQLite error";sqlite3_free(e);throw std::runtime_error(m);}}
int main(int argc,char** argv){
 sqlite3* db=nullptr;sqlite3_stmt* insert=nullptr;
 try{
  if(argc!=4)throw std::runtime_error("usage: gbz_graph_index graph.gbz NEW.sqlite metadata.json");
  std::ifstream existing(argv[2]);if(existing.good())throw std::runtime_error("Output already exists");
  std::ifstream meta_input(argv[3]);auto metadata=nlohmann::json::parse(meta_input);
  auto t=Clock::now();std::ifstream input(argv[1],std::ios::binary);
  if(!input)throw std::runtime_error("Cannot open GBZ");
  gbwtgraph::GBZ gbz;gbz.simple_sds_load(input);double load_seconds=elapsed(t);
  if(!gbz.index.bidirectional() || gbz.index.sequences()%2)throw std::runtime_error("Requires bidirectional GBWT");
  uint64_t paths=gbz.index.sequences()/2,max_id=gbz.graph.max_node_id(),nodes=gbz.graph.get_node_count();
  if(paths>INT32_MAX)throw std::runtime_error("Path counts may exceed tensor int32 capacity");
  bool dense=max_id<=4*nodes+1024 && max_id<=1000000000ULL;
  std::vector<uint32_t> counts,last;
  struct Entry{uint32_t count=0,last=0;};std::unordered_map<uint64_t,Entry> sparse;
  if(dense){counts.resize(max_id+1);last.resize(max_id+1);}else{sparse.reserve(nodes);}
  std::cerr<<"Loaded GBZ: "<<nodes<<" nodes, "<<paths<<" logical paths, "<<load_seconds<<" s\n";
  t=Clock::now();uint64_t visits=0;
  // Each logical path is walked once. last-seen tokens deduplicate revisits
  // including mixed orientations, without retaining entire extracted paths.
  for(uint64_t path=0;path<paths;++path){
   auto pos=gbz.index.start(gbwt::Path::encode(path,false));uint32_t token=path+1;
   while(pos.first!=gbwt::ENDMARKER){
    uint64_t id=gbwt::Node::id(pos.first);
    if(!id || id>max_id)throw std::runtime_error("GBWT node outside graph ID range");
    if(dense){if(last[id]!=token){last[id]=token;++counts[id];}}
    else{auto& e=sparse[id];if(e.last!=token){e.last=token;++e.count;}}
    ++visits;pos=gbz.index.LF(pos);
   }
   if((path+1)%100==0)std::cerr<<"Counted "<<path+1<<"/"<<paths<<" paths; "<<visits<<" visits; "<<elapsed(t)<<" s\n";
  }
  double count_seconds=elapsed(t);t=Clock::now();last.clear();last.shrink_to_fit();
  if(sqlite3_open_v2(argv[2],&db,SQLITE_OPEN_READWRITE|SQLITE_OPEN_CREATE,nullptr)!=SQLITE_OK)throw std::runtime_error("Cannot create SQLite");
  sql(db,"PRAGMA journal_mode=DELETE; PRAGMA synchronous=FULL; PRAGMA cache_size=-65536; CREATE TABLE graph_metadata(value TEXT NOT NULL); CREATE TABLE nodes(node_id INTEGER PRIMARY KEY CHECK(node_id>0), seq TEXT NOT NULL CHECK(length(seq)>0), distinct_path_count INTEGER NOT NULL CHECK(distinct_path_count>=0 AND distinct_path_count<=2147483647)); BEGIN;");
  if(sqlite3_prepare_v2(db,"INSERT INTO nodes VALUES(?,?,?)",-1,&insert,nullptr)!=SQLITE_OK)throw std::runtime_error(sqlite3_errmsg(db));
  uint64_t written=0;
  gbz.graph.for_each_handle([&](const handlegraph::handle_t& h){
   auto id=gbz.graph.get_id(h);auto seq=gbz.graph.get_sequence(gbz.graph.get_handle(id,false));uint32_t count=0;
   if(dense)count=counts.at(id);else{auto p=sparse.find(id);if(p!=sparse.end())count=p->second.count;}
   sqlite3_bind_int64(insert,1,id);sqlite3_bind_text(insert,2,seq.data(),seq.size(),SQLITE_TRANSIENT);sqlite3_bind_int64(insert,3,count);
   if(sqlite3_step(insert)!=SQLITE_DONE)throw std::runtime_error(sqlite3_errmsg(db));
   sqlite3_reset(insert);sqlite3_clear_bindings(insert);
   if(++written%100000==0){sql(db,"COMMIT; BEGIN;");if(written%1000000==0)std::cerr<<"Wrote "<<written<<" nodes\n";}
  });
  sqlite3_finalize(insert);insert=nullptr;if(written!=nodes)throw std::runtime_error("Graph node count mismatch");
  metadata["nodes"]=written;metadata["logical_paths"]=paths;metadata["path_visits"]=visits;metadata["status"]="complete";
  metadata["count_algorithm"]="one-forward-sequence-per-logical-path-last-seen-token-v1";
  metadata["timing"]={{"gbz_load_seconds",load_seconds},{"path_count_seconds",count_seconds},{"sqlite_write_seconds",elapsed(t)}};
  std::string value=metadata.dump();sqlite3_stmt* m=nullptr;
  if(sqlite3_prepare_v2(db,"INSERT INTO graph_metadata VALUES(?)",-1,&m,nullptr)!=SQLITE_OK)throw std::runtime_error(sqlite3_errmsg(db));
  sqlite3_bind_text(m,1,value.data(),value.size(),SQLITE_TRANSIENT);
  if(sqlite3_step(m)!=SQLITE_DONE){sqlite3_finalize(m);throw std::runtime_error(sqlite3_errmsg(db));}
  sqlite3_finalize(m);sql(db,"COMMIT;");sqlite3_close(db);db=nullptr;
  std::cout<<metadata.dump(2)<<std::endl;return 0;
 }catch(const std::exception& e){if(insert)sqlite3_finalize(insert);if(db)sqlite3_close(db);std::cerr<<"GBZ index build failed: "<<e.what()<<std::endl;return 1;}
}
