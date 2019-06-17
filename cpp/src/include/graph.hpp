#ifndef __BLOCK__GRAPH_HH
#define __BLOCK__GRAPH_HH


#include "intstack.hpp"

namespace block
{

template<class WEIGHT>
	struct weighted_edge {
		int vertex_;
		WEIGHT weight_;
		
		WEIGHT weight() { return weight_; }
		WEIGHT vertex() { return vertex_; }
		
		bool operator<(const weighted_edge& e) const 
		{
			return (weight_ < e.weight_) or (weight_ == e.weight_ and vertex_ > e.vertex_);
		}
		bool operator>=(const weighted_edge& e) const 
		{
			return !this->operator<(e);
		}
		bool operator==(const weighted_edge& e) const 
		{
			return weight_ == e.weight_ and vertex_ == e.vertex_;
		}
		bool operator>(const weighted_edge& e) const 
		{
			return e.operator<(*this);
		}
		bool operator<=(const weighted_edge& e) const 
		{
			return !e.operator>(*this);
		}
		bool operator!=(const weighted_edge& e) const 
		{
			return !this->operator==(e);
		}
		
	};


	template<class WEIGHT>
class graph
{

public:
    intstack nodes;
    std::vector<std::vector<weighted_edge<WEIGHT>>> neighbors;
		std::vector<WEIGHT> weight;

    unsigned int num_edges;

    graph() {}
    explicit graph(int n)
        : neighbors(n)
				, weight(n)
        , num_edges(0)
    {
				initialise(n);
		}
		
		void initialise(const int n) {
        nodes.reserve(n);
        nodes.fill();
    }

    size_t size() const { return nodes.size(); }
    size_t capacity() const { return nodes.capacity(); }

		void add_node(const int u, const WEIGHT w);
    void add_edge(const int u, const int v, const WEIGHT w);

    void sort();
		
		
		// void BFS(intstack& order);
		// template<class matrix, typename measure>
		// void FloydWarshall(matrix& distance, measure length);
		
		
    std::ostream& describe(std::ostream& os, const int verbosity) const;

    void check_consistency();
};


template<class WEIGHT>
	void graph<WEIGHT>::add_node(const int u, const WEIGHT w)
	{
	    weight[u] = w;
	}

template<class WEIGHT>
void graph<WEIGHT>::add_edge(const int u, const int v, const WEIGHT w)
{
    neighbors[u].push_back(weighted_edge<WEIGHT>{v,w});
    neighbors[v].push_back(weighted_edge<WEIGHT>{u,w});
    ++num_edges;
}

template<class WEIGHT>
void graph<WEIGHT>::sort()
{
    for (auto v : nodes)
        std::sort(begin(neighbors[v]), end(neighbors[v]), std::greater<weighted_edge<WEIGHT>>());
		std::vector<int> buffer;
		for (auto v : nodes)
				buffer.push_back(v);
		std::sort(begin(buffer), end(buffer), [&](int a, int b) { return weight[a] > weight[b]; });
		nodes.clear();
		for( auto v : buffer )
				nodes.add(v);
}

template<class WEIGHT>
std::ostream& graph<WEIGHT>::describe(std::ostream& os, const int verbosity) const
{
    // void gc::graph::describe(std::ostream& os, const int verbosity) const
    // {
    switch (verbosity) {
    case 0:
        os << "n=" << nodes.size() << " m=" << num_edges << std::endl;
        break;

    case 1:
        os << "V=" << nodes << " m=" << num_edges;
        break;

    case 2:
        for (auto v : nodes) {
            os << v << ":";
            for (auto u : neighbors[v]) {
                os << " " << u.vertex();
            }
            os << std::endl;
        }
  			break;
    case 3:
        for (auto v : nodes) {
            os << v << " (" << weight[v] << "):";
            for (auto u : neighbors[v]) {
                os << " " << u.vertex() << "/" << u.weight();
            }
            os << std::endl;
        }
    }

    return os;
}

template<class WEIGHT>
void graph<WEIGHT>::check_consistency()
{
    auto count{2 * num_edges};
    std::vector<int> f(nodes.capacity(), 0);
    std::vector<int> b(nodes.capacity(), 0);

    for (auto u : nodes) {
        for (auto v : neighbors[u]) {
            assert(nodes.contain(v));
            ++f[u];
            ++b[v];
            --count;
        }
    }

    if (count)
        std::cout << count << std::endl;

    assert(count == 0);

    for (auto u : nodes) {
        assert(f[u] == b[u]);
    }
}

template<class WEIGHT>
std::ostream& operator<<(std::ostream& os, const graph<WEIGHT>& g)
{
    return g.describe(os, 3);
}

} // namespace block

#endif
