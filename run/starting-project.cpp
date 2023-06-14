// Copyright © 2021 Giorgio Audrito. All Rights Reserved.


// [INTRODUCTION]
//! Importing the FCPP library.
#include "lib/fcpp.hpp"

const int NODE_NUM = 100;

namespace fcpp {

namespace coordination {

namespace tags {
    struct node_color {};
    struct node_size {};
    struct node_shape {};


    struct node_is_source_like{};
    struct node_is_sink_like{};

    struct flow_field {};
    struct node_excess {};

    //fields of outgoing capacities
    struct capacity_field {};
    struct residual_capacity_field {};


    struct node_distance_from_source{};
    struct node_distance_to_sink{};

}

//! @brief The maximum communication range between nodes.
constexpr size_t communication_range = 100;

const size_t dim = 2;



/*
In the description of the algorithm, weighted directed graphs on the set of devices will be 
identified with fields of reals:
if  d1 is a device, d2 a neighbour of d1 and the field F in d1  has value w nonzero in d2,
then we say that the graph associated with F has an arc from d1 to d2 with weight w.

This algorithm is inspired to the Ford–Fulkerson method: as long as there is an admissible 
path from source to sink in the residual graph, we push flow along the path.
We detect the presence of such paths simply by considering the distances from source and 
to sink along the residual graph.
If a node has both finite distance from source and finite distance to sink, then there is an
admissible path passing through it.
*/


/* All graphs considered here are directed. For this reason we have to implement a
 "distance along graph" function. We also are interested in which direction communicates to
  us which distance. Hence this function returns a field type.
*/
FUN field<real_t> distance_hood(ARGS, bool b, field<real_t>&& graph){ CODE

    return nbr(CALL, b?field<real_t>(0.0):field<real_t>(INF), [&](field<real_t> distances){
            real_t s = self(CALL, distances);
            

            field<real_t> tmp = map_hood([&](real_t d, real_t g){
                return g>0 ? d : INF;
            }, distances, graph);
            real_t m = min_hood(CALL, tmp);
    /*         
    our node communicates to neighbours its (guess regarding the) distance to nodes
    with a certain property. 
    In residual graphs there could be cycles that become disconnected from those nodes.
    Here there's a trick aimed at detecting those cycles. It doesn't work very well.
    I would need some more time in order to implement a message-efficient solution.
    */
            return make_tuple(distances, field<real_t>(b
                                                        ? 0.0
                                                        : s> m 
                                                            ? m+1
                                                            : INF));
    });
}

/* Here we implement our rule to assign dynamically the capacities between nodes.
These capacities will remain unaffected by the updates in our flow.
The Capacity weighted directed graph is acyclic.

Currently there are only two (similar) capacity functions proposed. This has capacities 
equal to 0 or 1.
*/
FUN field<real_t> capacity_v1(ARGS){ CODE
    field<device_t> ids = nbr_uid(CALL);
    return map_hood([&](device_t id){ return node.uid<id ;}, ids);
}

//This capacity funtion has different nonnull capacities.
FUN field<real_t> capacity_v2(ARGS){ CODE
    field<device_t> ids = nbr_uid(CALL);
    return map_hood([&](device_t id){ return node.uid<id ? id-node.uid:0;}, ids);
}

// Rough method to switch between capacity_v1 and capacity_v2
FUN field<real_t> capacity(ARGS){ CODE
    //return capacity_v1(CALL);
    return capacity_v2(CALL);
}


// Usual definition for the residual graph. Notice that it has only nonnegative weights.
// As already mentioned, it can have cycles.
FUN field<real_t> residual_capacity(ARGS, field<real_t> flow){ CODE
    return capacity(CALL) - flow;
}

//Reverse edges in a graph
FUN field<real_t> incoming(ARGS, field<real_t>&& graph){ CODE
    return nbr(CALL, graph);
}

// Reverse the residual graph
FUN field<real_t> incoming_residual_capacity(ARGS, field<real_t> flow){ CODE
    return incoming(CALL, residual_capacity(CALL, flow));
}

//This is just another function to sum up the values in a field
real_t sum(field<real_t>& input){

    real_t tmp= 0.0;
    for (details::field_iterator<field<real_t>> it(input); !it.end(); ++it){
        tmp += it.value();
    }
    return tmp;
}

/*
Function that monitors the ability of a node to push or receive flow.
Strictly negative excess means that a node receives more flow in ingress so it can push some.
Source is always able to push flow.
Strictly positive excess means that a node pushes more flow so it has to receive some.
Sink is always ready to receive flow.
*/
FUN real_t excess(ARGS, field<real_t> flow){ CODE

    bool is_source = node.uid==0;
    bool is_sink = node.uid == NODE_NUM-1;
    return is_source
                ?-INF
                :is_sink
                    ? INF
                    :sum( flow);
}


// We call source-like a node that has to push flow
FUN bool is_source_like(ARGS, field<real_t> flow){ CODE
    return excess(CALL, flow)<0;
}

// We call sink-like a node that has to receive flow
FUN bool is_sink_like(ARGS, field<real_t> flow){ CODE
    return excess(CALL, flow)>0;
}

// Returns the field of distances from the closest source-like nodes
FUN field<real_t> from_source_field(ARGS, field<real_t> flow){ CODE
    return distance_hood(CALL, is_source_like(CALL, flow), incoming_residual_capacity(CALL, flow));
}

// Returns the distance from the closest source-like node
FUN real_t from_source(ARGS, field<real_t> flow){ CODE
    return self(CALL, from_source_field(CALL, flow));
}

// Returns the field of distances to the closest sink-like nodes
FUN field<real_t> to_sink_field(ARGS, field<real_t> flow){ CODE
    return distance_hood(CALL, is_sink_like(CALL, flow), residual_capacity(CALL, flow));
}

// Returns the distance to the closest sink-like node
FUN real_t to_sink(ARGS, field<real_t> flow){ CODE
    return self(CALL, to_sink_field(CALL, flow));
}

// Returns the increment of flow per round
FUN field<real_t> flow_increment(ARGS, field<real_t> flow){ CODE
    // _n stands for "name"
    
    field<real_t> residual_capacity_n = residual_capacity(CALL, flow);

    field<real_t> to_sink_field_n = to_sink_field(CALL, flow);
    real_t to_sink_n = to_sink(CALL, flow);

    real_t excess_n = excess(CALL, flow);

    // In this case we push flow along minimal admissible paths, 
    //until we have flow in ingress to push forward
    if(to_sink_n!= INF  && excess_n<0){
        return map_hood([&](real_t d, real_t r){
                                        real_t a = 0.0;
                                        if(excess_n<0 && d == to_sink_n-1){
                                            a = std::min(r, -excess_n);
                                            excess_n+= a;
                                        }
                                        return a;
                                    },to_sink_field_n, residual_capacity_n);
    }
    // In this case the node realizes to not be anymore in an admissible path 
    // and so we have to push back some flow
    else if(excess_n<0){
        return map_hood([&](real_t f){
                        real_t a = 0.0;
                        if(excess_n<0 && f<0){
                            a = std::min(-f, -excess_n);
                            excess_n += a;
                        }
                        return a;
                    }, flow); 
    }
    // Nodes that don't have flow to push, don't increment flow.
    else return field<real_t>(0.0);
                
}


//Updates the flow adding the increment
FUN field<real_t> update_flow(ARGS){ CODE
    return nbr(CALL, field<real_t>(0.0),[&](field<real_t> flow){
        field<real_t> up_flow = flow - nbr(CALL, flow_increment(CALL, flow));
        return make_tuple(up_flow, -up_flow);
    });
}







//! @brief Main function.
MAIN() {
    // import tag names in the local scope.
    using namespace tags;

    // References
    field<real_t>& capacity_ = node.storage(capacity_field{});
    field<real_t>& flow_ = node.storage(flow_field{});
    field<real_t>& residual_capacity_ = node.storage(residual_capacity_field{});
    real_t& from_source_ = node.storage(node_distance_from_source{});
    real_t& to_sink_ = node.storage(node_distance_to_sink{});
    real_t& excess_ = node.storage(node_excess{});
    bool& is_source_like_ = node.storage(node_is_source_like{});
    bool& is_sink_like_ = node.storage(node_is_sink_like{});


    // Usage of node storage
    node.storage(node_size{}) = 10;
   

    bool is_source = node.uid==0;
    bool is_sink = node.uid == NODE_NUM-1;

    node.storage(node_shape{}) = is_source
                                    ?shape::star 
                                    :is_sink
                                        ? shape::tetrahedron
                                        : shape::sphere;

    // This is the only structure that node requires to manage
    flow_ = update_flow(CALL);
    
    
    // These other structures are just aimed at monitoring the behaviour of system
    capacity_ = capacity(CALL);
    residual_capacity_ = residual_capacity(CALL, flow_);
    is_source_like_ = is_source_like(CALL, flow_);
    is_sink_like_ = is_sink_like(CALL, flow_);
    from_source_ = from_source(CALL, flow_);
    to_sink_ = to_sink(CALL, flow_);
    /*
    In this structurre we monitor how much flow source pushes and
    how much flow sink receives. Hopefully they're equal in absolute module
    */
    excess_= sum( flow_);

    /*
    Nodes that have flow to push are GREEN; those that need to receive flow are RED;
    Nodes that guess to be part of an admissible path from a source-like to a sink-like
    are YELLOW; other nodes are WHITE.
    */
    node.storage(node_color{}) =   is_source_like_
                                    ? color(GREEN)
                                    : is_sink_like_
                                        ? color(RED)
                                        : from_source_!=INF && to_sink_!=INF
                                            ?color(YELLOW)
                                            :color(WHITE);
}
//! @brief Export types used by the main function.
FUN_EXPORT main_t = export_list<device_t, field<real_t>, real_t, bool >;

} // namespace coordination

// [SYSTEM SETUP]

//! @brief Namespace for component options.
namespace option {

//! @brief Import tags to be used for component options.
using namespace component::tags;
//! @brief Import tags used by aggregate functions.
using namespace coordination::tags;

//! @brief Number of people in the area.
constexpr int node_num = NODE_NUM;
//! @brief Dimensionality of the space.
constexpr size_t dim = 2;

//! @brief Description of the round schedule.
using round_s = sequence::periodic<
    distribution::interval_n<times_t, 0, 1>,    // uniform time in the [0,1] interval for start
    distribution::weibull_n<times_t, 10, 1, 10> // weibull-distributed time for interval (10/10=1 mean, 1/10=0.1 deviation)
>;
//! @brief The sequence of network snapshots (one every simulated second).
using log_s = sequence::periodic_n<1, 0, 1>;
//! @brief The sequence of node generation events (node_num devices all generated at time 0).
using spawn_s = sequence::multiple_n<node_num, 0>;
//! @brief The distribution of initial node positions (random in a 500x500 square).
using rectangle_d = distribution::rect_n<1, 0, 0, 500, 500>;
//! @brief The contents of the node storage as tags and associated types.
using store_t = tuple_store<
    node_color,                         color,
    node_size,                          double,
    node_shape,                         shape,
    node_is_source_like,                bool,
    node_is_sink_like,                  bool,
    node_distance_from_source,          real_t,
    node_distance_to_sink,              real_t,
    node_excess,                        real_t,
    capacity_field,                     field<real_t>,
    flow_field,                         field<real_t>,
    residual_capacity_field,            field<real_t>

>;
//! @brief The tags and corresponding aggregators to be logged (change as needed).
using aggregator_t = aggregators<
>;

//! @brief The general simulation options.
DECLARE_OPTIONS(list,
    parallel<true>,      // multithreading enabled on node rounds
    synchronised<false>, // optimise for asynchronous networks
    program<coordination::main>,   // program to be run (refers to MAIN above)
    exports<coordination::main_t>, // export type list (types used in messages)
    retain<metric::retain<2,1>>,   // messages are kept for 2 seconds before expiring
    round_schedule<round_s>, // the sequence generator for round events on nodes
    log_schedule<log_s>,     // the sequence generator for log events on the network
    spawn_schedule<spawn_s>, // the sequence generator of node creation events on the network
    store_t,       // the contents of the node storage
    aggregator_t,  // the tags and corresponding aggregators to be logged
    init<
        x,      rectangle_d // initialise position randomly in a rectangle for new nodes
    >,
    dimension<dim>, // dimensionality of the space
    connector<connect::fixed<100, 1, dim>>, // connection allowed within a fixed comm range
    shape_tag<node_shape>, // the shape of a node is read from this tag in the store
    size_tag<node_size>,   // the size  of a node is read from this tag in the store
    color_tag<node_color>  // the color of a node is read from this tag in the store
);

} // namespace option

} // namespace fcpp


//! @brief The main function.
int main() {
    using namespace fcpp;

    //! @brief The network object type (interactive simulator with given options).
    using net_t = component::interactive_simulator<option::list>::net;
    //! @brief The initialisation values (simulation name).
    auto init_v = common::make_tagged_tuple<option::name>("Exercises");
    //! @brief Construct the network object.
    net_t network{init_v};
    //! @brief Run the simulation until exit.
    network.run();
    return 0;
}
