window.FabricConfig = { fontBaseUrl: 'https://your-url-here' }

function showtree(drawing_type, dst_id, newick, width, height)
{
	console.log("showtree :: drawing_type: " + drawing_type + " dst_id: " + dst_id + " width: " + width + " height: " + height);

    var radial = drawing_type == "circle";

    // http://phylotree.hyphy.org/documentation/options.html#options-section
    var tree = d3.layout.phylotree()
        .svg(d3.select("#"+dst_id))
        .options({
            'collapsible': true,
            'zoom': true,
            'align-tips': true,
            'brush': false
        })
        .size([height, width])
        .node_circle_size(0)
        .radial(radial);
        ;
    
    tree(newick)
        .size([height, width])
        .font_size(12)
        .layout();
}