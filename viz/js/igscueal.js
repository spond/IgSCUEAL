var master_array;
var chart_dim = 400;
var jd_handler;
var is_light_chain = false;
var j_column_id = 3;

function objectType (obj) {
    return Object.prototype.toString.call(obj).split(' ').pop().split(']').shift().toLowerCase();
}

function callJD (event) {
    return jd_handler (event);
}

function jdClicker (event, label) {
    ed = d3.select (event.srcElement).data();
    draw_svg_barchart("cdr3_chart", bin_by_cdr3 (master_array,function (d) {return d[1] == label && d[2] == ed[0][0] && d[3] == ed[0][1];}), 300, 300, label+","+ed[0][0]+","+ed[0][1]+ " rearrangement", "Length", "Count","cdr3_table");
}

function tertiaryChart (event) {
    ed = d3.select (event.srcElement).data();
    if (objectType (ed[0]) != 'array') {
        var filter_on = ed[0];
    } else {
        var filter_on = ed[0][0];
    }
    var filter = function (d) {return d[1] == filter_on;};
    
    jd_handler = function  (event) {
        jdClicker (event, filter_on);
    }
                    
    draw_svg_scatter("jdchart", bin_array_by_v_j (master_array,2,3,filter),chart_dim,chart_dim, filter_on + " allele", "D-allele", "J-allele", [callJD,null,null]);
    draw_svg_barchart("cdr3_chart", bin_by_cdr3 (master_array,filter), 300, 300, filter_on + " allele", "Length", "Count","cdr3_table");
}

function secondaryChart (event) {
    ed = d3.select (event.srcElement).data();
    if (objectType (ed[0]) != 'array') {
        var filter_on = ed[0];
    } else {
        var filter_on = ed[0][0];
    }
    var filter = function (d) {return d[0] == filter_on;};        

    draw_svg_scatter("v2jchart", bin_array_by_v_j (master_array,1,j_column_id,filter),chart_dim,chart_dim, filter_on + " family", "V-allele", "J-allele", is_light_chain?[null,null,null]:[tertiaryChart,tertiaryChart,null]);
    draw_svg_barchart("cdr3_chart", bin_by_cdr3 (master_array,filter), 300, 300, filter_on + " family", "Length", "Count","cdr3_table");
}

function secondaryChartY (event) {
    if (!is_light_chain) {
        ed = d3.select (event.srcElement).data();
        var filter = function (d) {return d[3] == ed[0];};
        draw_svg_scatter("jdchart", bin_array_by_v_j (master_array,1,2,filter),chart_dim,chart_dim, ed[0] + " allele", "V-allele", "D-allele", "");
        draw_svg_barchart("cdr3_chart", bin_by_cdr3 (master_array,filter), 300, 300, ed[0] + " allele", "Length", "Count","cdr3_table");
    }
}

function hist_mean (array, total) {
    var mean = 0.;
    for (k in array) {
        mean += array[k][0] * array[k][1]/total;
    }   
    return mean;
}

function hist_quantile (array, total, p) {
    if (p == 1) {
        return array[array.length-1][0];
    }
    if (array.length == 1) {
        return array[0][0];
    }
    
    var h    = (total-1)*p;
    var hint = Math.floor (h);
    var i = 0;
    
    var sum = array[0][1]-1;
    while (sum < hint) {
        i++;
        sum += array[i][1];
    }
    
    var i2 = i;
    while (sum < h) {
        i2++;
        sum += array[i2][1];
    }
            
    if (i == i2) {
        return array[i][0];
    } 
    
    return array[i][0] + (h-hint)*(array[i2][0]-array[i][0]);
}

function make_stats (array) {
    array.sort(function (row1, row2) { return row1[0] - row2[0];});
    var result = [];
    var total  = d3.sum (array, function (d) {return d[1];});
    result.push (['Total', total]);
    result.push (['Mean',  hist_mean(array,total).toFixed(2)]);
    result.push (['Median',  hist_quantile (array, total, 0.5).toFixed(2)]);
    result.push (['Min',  array[0][0].toFixed(2)]);
    result.push (['1.0%',  hist_quantile (array, total, 0.01).toFixed(2)]);
    result.push (['2.5%',  hist_quantile (array, total, 0.025).toFixed(2)]);
    result.push (['97.5%',  hist_quantile (array, total, 0.975).toFixed(2)]);
    result.push (['99.0%',  hist_quantile (array, total, 0.99).toFixed(2)]);
    result.push (['Max',  array[array.length-1][0].toFixed(2)]);
    
    return result;
}

function make_frequency_data (array, column) {
    var result = {};
    array.sort(function (row1, row2) { if(row1[column] < row2[column]) return -1;if(row1[column] > row2[column]) return 1;return 0;});
    var return_me = [['Assignment', 'Count']];

    for (k in array) {
        var key = array[k][column];
        if (! (key in result)) {
            result [key] = array[k][2];
        } else {
            result [key] += array[k][2];
        }
    }   
    for (var key in result) {
        return_me.push([key, result[key]]);
    }
    return return_me;
}

function draw_svg_barchart (element_tag, input_array, width,height, chartLabel, labelX, labelY, tableTag, doLog) {
    var data = input_array["Lengths"];
    var x_scale, y_scale;
    if (doLog) {
        x_scale = d3.scale.log().domain ([Math.min.apply (0, data.map (function (d) {return d[0]})),Math.max.apply (0, data.map (function (d) {return d[0];}))]);
        y_scale = d3.scale.log().domain ([Math.max.apply (0, data.map (function (d) {return d[1];})),0]);
    
    } else {
        x_scale = d3.scale.linear().domain ([Math.min.apply (0, data.map (function (d) {return d[0]})),Math.max.apply (0, data.map (function (d) {return d[0];}))]);
        y_scale = d3.scale.linear().domain ([Math.max.apply (0, data.map (function (d) {return d[1];})),0]);
    }
    
    var max_count = y_scale.domain()[1];


    var chart_padding = 40;
    var title_padding = 20;

    x_scale.range ([chart_padding,width-chart_padding]);
    y_scale.range ([title_padding,height-chart_padding]);
    var bar_width = x_scale(max_count) - x_scale (max_count - 1);

    var chart = d3.select("#"+element_tag)
                .attr ("width",  width)
                .attr ("height", height);
                

    chart.selectAll ("rect").data([]).exit().remove();
    chart.selectAll ("g").remove();

    var bars = chart.selectAll ("rect").data (data);
    var y_span = y_scale.range();
    y_span = y_span[1] - y_span[0];

    bars  = bars.enter()
                 .append ("rect")
                 .attr ('x', function (d) {return x_scale(d[0]);})
                 .attr ('y', function (d) {return y_scale(d[1]);})
                 .attr ('width',  bar_width)
                 .attr ('height', function (d) {return y_span+title_padding-y_scale(d[1]);})
                 .attr ('class', 'cdr3bar')
                 .append ('title')
                 .text (function (d) {return d[0] + " (" + d[1] + ")";});

    var x_axis = d3.svg.axis().scale(x_scale).orient ("bottom");
    var y_axis = d3.svg.axis().scale(y_scale).orient ("left");
    var x_g = chart.append("g")
          .attr("class", "axis")
          .attr("transform", "translate("+bar_width/2+"," + (height-chart_padding+1) + ")")
          .call(x_axis);   
             
     x_g = chart.append("g")
          .attr("class", "axis")
          .attr("transform", "translate(" + chart_padding +",0)")
          .call(y_axis) ;      

     x_g.selectAll ("text").attr ("transform", "rotate(-45)").attr("dx","0.5em").attr("dy","-1em");
     var label = chart.selectAll (".charttitle");
     if (label.empty()) {
        label = chart.append ("text").attr ("class", "charttitle");
     }
     
          
     label.attr ("text-anchor", "middle")
            .attr ("x",width/2)
            .attr ("y",title_padding-5)
            .text (chartLabel);     
            
     if (tableTag.length > 0) {
        plot_table_with_frequencies (tableTag,make_stats(data), false, "table_cell_small");
    }

}

function draw_svg_scatter (element_tag, paired_array, width,height, chartLabel, labelX, labelY, clicker) {
        
    var x_scale = d3.scale.ordinal().domain (paired_array.map (function (d) {return d[0]}).sort());
    var y_scale = d3.scale.ordinal().domain (paired_array.map (function (d) {return d[1]}).sort());
    
    var max_count = Math.max.apply (0, paired_array.map (function (d) {return d[2]}));
    var ref_label = Math.pow(10,Math.floor(Math.log(max_count)/Math.log(10)));
    var reference_size = ref_label/max_count;
    var max_v     = Math.sqrt(Math.max.apply (0, paired_array.map (function (d) {return d[3]})));
    
    var total     = d3.sum (paired_array, function (d) {return d[2];});

    var chart_padding = 35;
    var title_padding = 20;
    var legend_width  = 100;
    
    width = Math.max (width, 15*x_scale.domain ().length);
    x_scale.rangeRoundBands ([chart_padding*2,width], 0.0, 0.5);
    y_scale.rangeRoundBands ([chart_padding+title_padding,height-chart_padding],0.0, 0.5);
    
    
    var x_band = x_scale.rangeBand();
    var y_band = y_scale.rangeBand();
    
    radius = Math.min (x_band,y_band)/2;
                    
    var chart = d3.select("#"+element_tag)
                .attr ("width",  width)
                .attr ("height", height);
     
    var x_axis = d3.svg.axis().scale(x_scale).orient ("bottom");
    var y_axis = d3.svg.axis().scale(y_scale).orient ("left");
             
    x_axis (chart);
    y_axis (chart);
    
    chart.selectAll ("circle").data([]).exit().remove();
    chart.selectAll ("g").remove();
                              
    var circles = chart.selectAll ("circle").data (paired_array);
    

    circles = circles.enter()
                 .append ("circle")
                 .attr ('cx', function (d) {return x_scale(d[0]);})
                 .attr ('cy', function (d) {return y_scale(d[1]);})
                 .attr ('r',  function (d) {return Math.max(1,radius*(Math.sqrt(d[3])/max_v));})
                 .attr ('class', 'igblob')
                 .append ('title')
                 .text (function (d) {return d[2] + " (" + (d[3]*100).toFixed(2) + "%)";});
                 
     if (clicker.length) {
        if (clicker[0] != null) {
            chart.selectAll ("circle").attr ('onclick', clicker[0].name + "(event)");
        }
     }            
          
                              
     var x_g = chart.append("g")
          .attr("class", "axis")
          .attr("transform", "translate(-"+x_band*0.5+"," + (height-chart_padding) + ")")
          .call(x_axis);      
          
     if (clicker.length > 1 && clicker[1] != null) {     
         x_g.selectAll("text").attr ("transform", "rotate(-45)").attr("dx","-1em").attr ('onclick', clicker[1].name + "(event)");
     } else {
         x_g.selectAll("text").attr ("transform", "rotate(-45)").attr("dx","-1em");
      }
     x_g = chart.append("g")
          .attr("class", "axis")
          .attr("transform", "translate(" + chart_padding +",-"+y_band*0.5+")")
          .call(y_axis) ;      
      
     if (clicker.length > 2 && clicker[2] != null) {     
        x_g.selectAll ("text").attr ("transform", "rotate(-45)").attr("dx","0.5em").attr("dy","-1em").attr ('onclick', clicker[2].name + "(event)");
     } else {
        x_g.selectAll ("text").attr ("transform", "rotate(-45)").attr("dx","0.5em").attr("dy","-1em");
     }

     var legend_blob = chart.append ("g").attr ("transform", "translate (" + (width - legend_width) + "," + 0 + ")");     
     var legend_blob_radius = radius * reference_size;
     legend_blob.append ("circle").attr ("cx", radius).attr("cy", radius).attr("r", legend_blob_radius).attr('class','iglegend');
     
     if (legend_blob_radius > 10) {
        legend_blob.append ("text").attr("text-anchor", "middle").attr("dominant-baseline", "central").attr ("x",radius).attr("y",radius).text (ref_label).attr('class','legend');            
     } else {
        legend_blob.append ("text").attr("text-anchor", "start").attr("dominant-baseline", "central").attr ("x",radius + legend_blob_radius + 5).attr("y",radius).text (ref_label).attr('class','legend');
     }
     
     var label = chart.selectAll (".charttitle");
     if (label.empty()) {
        label = chart.append ("text").attr ("class", "charttitle");
     }
     
          
     label.attr ("text-anchor", "middle")
            .attr ("x",width/2)
            .attr ("y",title_padding-5)
            .text (chartLabel);     

      label = chart.selectAll (".xaxislabel");
      if (label.empty()) {
          label = chart.append ("text").attr ("class", "xaxislabel");
      }
     
      label.attr("text-anchor", "end")
            .attr("x", width)
            .attr("y", height - chart_padding)
            .attr("dy", "-.75em")
            .text(labelX);
            
      label = chart.selectAll (".yaxislabel");
      if (label.empty()) {
         label = chart.append ("text").attr ("class", "yaxislabel");
      }

      label.attr("text-anchor", "end")
            .attr("y", 6)
            .attr("dy", ".75em")
            .attr("transform", "rotate(-90)")
            .text(labelY);
        
       plot_table_with_frequencies ("freqs_table_x",make_frequency_data (paired_array, 0), total, "table_cell_small");
       plot_table_with_frequencies ("freqs_table_y",make_frequency_data (paired_array, 1), total, "table_cell_small");
}


function bin_array_by_v_j (array, col1, col2, apply_filter) {
    var resulting_array = {}; // v, j, count
    var total_count = 0;
    for (key in array) {
        if (key != 'Summary') {
            var bits = key.split('|'); 
            if (apply_filter) {
                if (apply_filter (bits) == false) {
                    continue;
                }
            }
            total_count += array[key]['Count']; 
            bits = bits[col1] + '|'+ bits[col2]; 
             if (!(bits in resulting_array)) {
                resulting_array[bits] = array[key]['Count'];
            } else {
                resulting_array[bits] += array[key]['Count'];
            }
        }
    }
    
    binned = []
    for (key in resulting_array) {
        components = key.split ('|');
        binned.push ([components[0], components[1], resulting_array[key], resulting_array[key]/total_count])
    }
    return binned;
}

function bin_by_cdr3 (array, apply_filter) {
    var binning = {}; // length: count, 
    var total_count = 0;
    resulting_array = {"Lengths": []};
    for (key in array) {
        if (key != 'Summary') {
            var bits = key.split('|'); 
            if (apply_filter) {
                if (apply_filter (bits) == false) {
                    continue;
                }
            }
            total_count += array[key]['Count']; 
            
            for (value in array[key]) {
                if (value != 'Count') {
                    var cdr3l = value-2;
                    if (!(cdr3l in binning)) {
                        binning[cdr3l] = array[key][value];
                    } else {
                       binning[cdr3l] += array[key][value];                        
                    }
                }
            }                 
        }
    }
    
    for (key in binning) {
        resulting_array["Lengths"].push ([parseInt(key), binning[key]]);
    }
    
    resulting_array["Count"] = total_count;
    return resulting_array;
}

function plot_table_with_frequencies (tag, table_data, do_proportions, table_class) {
    var do_tab = false;
    var max_value = 0;
    if (objectType (do_proportions) == 'number') {
        max_value = do_proportions;
        do_tab = true;
    } else {
        if (do_proportions) {
            max_value = d3.max (table_data, function (row, index) { if (index == 0) {return 0.0;} return row[1];});
            do_tab = true;
        } 
    }
    if (do_tab) {          
        table_data = table_data.map (function (row, index) { if (index>0) {return [row[0], row[1], (row[1]/max_value*100.0).toFixed(2) + "%"];} else {return [row[0],row[1],''];} }) ;
        table_data[0].splice(-1,0,'Proportion');
    }
    
    d3.select ('#'+tag).style("display","inline").selectAll('tr').data([]).exit().remove();
    var all_rows = d3.select ('#'+tag).selectAll('tr').data(table_data).enter().append ("tr");
    all_rows.selectAll ("td").data (function (a_row) {return a_row;}).enter().append ("td").text (function (text_value) {return text_value}).attr ("class", table_class).
        style("background-color", function (d,i,j) {if (j%2 == 1) return "#cccccc"; return "white"});
}

function load_igg_result  (json_object) {
     summary = json_object["Summary"];
     summary.splice (0,0,['','Count']);
     master_array = json_object;
     draw_svg_scatter("vj_chart", bin_array_by_v_j (json_object,0,j_column_id),chart_dim,chart_dim, "All rearrangements", "V-family", "J-allele",[secondaryChart,secondaryChart,secondaryChartY]);
      
     draw_svg_barchart("cdr3_chart", bin_by_cdr3 (json_object,""), 300, 300, "All rearrangements", "Length", "Count","cdr3_table");
      
     //draw_svg_scatter(bin_array_by_v_j (json_object,1,3, function (d) {return d[0] == 'V1';}),400,400);
     plot_table_with_frequencies ('summary_table', summary, true, "table_cell");
     d3.select("#b_cell_display").style ("display", "inline");
     return 0;
}

function hide_charts (table_list) {
    for (var k in table_list) {
        var table_id = table_list[k];
        d3.select("#" + table_id).text("");
    }
}

function hide_tables (table_list) {
    for (var k in table_list) {
        var table_id = table_list[k];
        d3.select("#" + table_id).selectAll("tr").data([]).exit().remove();
    }
}
    
    
function load_result_from_choice () {
    var selector =  d3.select ('#result_selector').node();
    hide_charts (["vj_chart","v2jchart","jdchart","cdr3_chart"]);
    hide_tables (["summary_table","cdr3_stats","freqs_table_x","freqs_table_y"]);
    var run_info = d3.select (selector.options[selector.selectedIndex]).datum();
    is_light_chain = (run_info[2] == "light");
    j_column_id = 3-is_light_chain;
    d3.json(run_info[1], load_igg_result);
}

function load_result_options (json_object) {
    var selector = d3.select ('#result_selector'); 
    selector.selectAll ("option").data([]).exit().remove();
    selector.selectAll ("option").data(json_object).enter().append ("option")
                .text (function (d) {return d[0];});
    selector.node().onchange();      
}

function load_file (url, chain_type) {

    hide_charts (["vj_chart","v2jchart","jdchart","cdr3_chart"]);
    hide_tables (["summary_table","cdr3_stats","freqs_table_x","freqs_table_y"]);

    d3.json (url, function (data) {
        is_light_chain = chain_type == 'light' ? true : false;
        j_column_id = 3-is_light_chain;
        load_igg_result (data);
    });
}