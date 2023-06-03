// Copyright © 2021 Giorgio Audrito. All Rights Reserved.


// [INTRODUCTION]
//! Importing the FCPP library.
#include "lib/fcpp.hpp"

const int NODE_NUM = 100;

/**
 * @brief Namespace containing all the objects in the FCPP library.
 */
namespace fcpp {

//! @brief Dummy ordering between positions (allows positions to be used as secondary keys in ordered tuples).
template <size_t n>
bool operator<(vec<n> const& x, vec<n> const& y) {
    return x[0]<y[0];
}

//! @brief Namespace containing the libraries of coordination routines.
namespace coordination {

//! @brief Tags used in the node storage.
namespace tags {
    //! @brief Color of the current node.
    struct node_color {};
    //! @brief Size of the current node.
    struct node_size {};
    //! @brief Shape of the current node.
    struct node_shape {};
    // ... add more as needed, here and in the tuple_store<...> option below


    struct node_is_source_like{};
    struct node_is_sink_like{};

    struct flow_field {};
    struct node_excess {};
    //fields of outgoing capacities
    struct capacity_field {};
    struct residual_capacity_field {};

    //fields of incoming capacities
    struct incoming_residual_capacity_field {};

    //fields of distances
    struct distance_from_source_field{};
    //struct distance_to_sink_field{};

    struct node_distance_from_source{};
    struct node_distance_to_sink{};

}

//! @brief The maximum communication range between nodes.
constexpr size_t communication_range = 100;



const size_t dim = 2;




real_t sum(field<real_t>& input){

    real_t tmp= 0.0;
    for (details::field_iterator<field<real_t>> it(input); !it.end(); ++it){
        tmp += it.value();
    }
    return tmp;
}


real_t min(field<real_t> input){

    real_t tmp= INF;
    for (details::field_iterator<field<real_t>> it(input); !it.end(); ++it){
        tmp = std::min(it.value(), tmp);
    }
    return tmp;
}

FUN field<real_t> distance_hood(ARGS, bool b, field<real_t>& graph){ CODE

    return nbr(CALL, b?field<real_t>(0.0):field<real_t>(INF), [&](field<real_t> distances){
            real_t& self = mod_self(CALL, distances);
            
            field<real_t> tmp = map_hood([&](real_t d, real_t g){
                return g>0 ? d : INF;
            }, distances, graph);
            real_t m = min(tmp);

// funziona abbastanza, ma quando la rete è intricata, il flusso rimane instabile vicino alla sorgente
            self = b
                ? 0.0
                : self> m ? m+1: INF;

       
             
                                        
            return make_tuple(distances, field<real_t>(self));


    });
}









//! @brief Main function.
MAIN() {
    // import tag names in the local scope.
    using namespace tags;

    //references
    field<real_t>& capacity = node.storage(capacity_field{});
    field<real_t>& flow_ = node.storage(flow_field{});
    field<real_t>& residual_capacity = node.storage(residual_capacity_field{});
    field<real_t>& incoming_residual_capacity = node.storage(incoming_residual_capacity_field{});
    real_t& from_source = node.storage(node_distance_from_source{});
    real_t& to_sink = node.storage(node_distance_to_sink{});
    real_t& excess_ = node.storage(node_excess{});
    field<real_t>& from_source_field = node.storage(distance_from_source_field{});
    //field<real_t>& to_sink_field = node.storage(distance_to_sink_field{});
    bool& is_source_like = node.storage(node_is_source_like{});
    bool& is_sink_like = node.storage(node_is_sink_like{});


        // usage of node storage
    node.storage(node_size{}) = 10;
   
    bool is_source = node.uid==0;
    bool is_sink = node.uid == NODE_NUM-1;

    node.storage(node_shape{}) = is_source
                                    ?shape::star 
                                    :is_sink
                                        ? shape::tetrahedron
                                        : shape::sphere;

    field<device_t> ids = nbr_uid(CALL);
    //we want  that, for every couple of neighbours, capacity is nonzero for at most one direction
    // ? id-node.uid:0
    capacity = map_hood([&](device_t id){ return node.uid<id ? id-node.uid:0;}, ids); 

    flow_ = nbr(CALL, field<real_t>(0.0),[&](field<real_t> flow){


            real_t excess = is_source
                            ?-INF
                            :is_sink
                                ? INF
                                :sum(flow);

            is_source_like = excess<0;
            is_sink_like = excess>0;

            residual_capacity = capacity - flow;
            incoming_residual_capacity = nbr(CALL, residual_capacity);

            from_source_field = distance_hood(CALL, is_source_like, incoming_residual_capacity);
            from_source = self(CALL, from_source_field);


            to_sink = nbr(CALL, is_sink_like?0.0:INF , [&](field<real_t> distances){

                real_t& self = mod_self(CALL, distances);

                field<real_t> tmp = map_hood([&](real_t d, real_t r, real_t i){
                    return  i>0 && d<0
                                ? -1
                                : r>0 
                                    ? d 
                                    : INF;
                }, distances, residual_capacity, incoming_residual_capacity);
                real_t m = min(tmp);
                self = is_sink_like 
                    ? 0.0
                    : m<0
                        ? INF
                        : self>m
                            ? m+1
                            : INF;

                
                field<real_t> flow_increment(0.0);
                real_t exc = excess;

                if(self!= INF  && exc<0){
                    flow_increment = map_hood([&](real_t d, real_t r){
                                        real_t a = 0.0;
                                        if(exc<0 && d == self-1){
                                            a = std::min(r, -exc);
                                            exc+= a;
                                        }
                                        return a;
                                    },distances, residual_capacity);
                
                }
                if(self==INF && exc<0){
                    flow_increment += map_hood([&](real_t f){
                        real_t a = 0.0;
                        if(exc<0 && f<0){
                            a = std::min(-f, -exc);
                            exc += a;
                        }
                        return a;
                    }, flow); 
                }

                flow -= nbr(CALL, flow_increment);
                
                return self;
            });

            
        return make_tuple(flow, -flow);
    });
    
    excess_= sum(flow_);


    node.storage(node_color{}) =   is_source_like
                                    ? color(GREEN)
                                    : is_sink_like
                                        ? color(RED)
                                        : from_source!=INF && to_sink!=INF
                                            ?color(YELLOW)
                                            :color(WHITE);

   //node.velocity() = -node.position()/300;
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
    residual_capacity_field,            field<real_t>,
    incoming_residual_capacity_field,   field<real_t>,
    distance_from_source_field,         field<real_t>
    //distance_to_sink_field,             field<real_t>

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
