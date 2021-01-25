// http://stackoverflow.com/questions/498970/how-do-i-trim-a-string-in-javascript
if (!String.prototype.trim)
{
	String.prototype.trim=function(){return this.replace(/^\s+|\s+$/g, '');};
}

function showtree(drawing_type, dst_id, newick, width, height)
{
	var t   = new Tree();
	newick  = newick.trim(newick);
	t.Parse(newick);

	if (t.error != 0)
	{
		document.getElementById('message').innerHTML='Error parsing tree';
	}
	else
	{
		document.getElementById('message').innerHTML='Parsed OK';
		
		t.WriteNewick();
		
		t.ComputeWeights(t.root);
		
		var td = null;
		
		// var selectmenu = document.getElementById('style');
		// var drawing_type = (selectmenu.options[selectmenu.selectedIndex].value);
		
		switch (drawing_type)
		{
			case 'rectanglecladogram':
				td = new RectangleTreeDrawer();
				break;
		
			case 'phylogram':
				if (t.has_edge_lengths)
				{
					td = new PhylogramTreeDrawer();
				}
				else
				{
					td = new RectangleTreeDrawer();
				}
				break;
				
			case 'circle':
				td = new CircleTreeDrawer();
				break;
				
			case 'circlephylogram':
				if (t.has_edge_lengths)
				{
					td = new CirclePhylogramDrawer();
				}
				else
				{
					td = new CircleTreeDrawer();
				}
				break;
				
			case 'cladogram':
			default:
				td = new TreeDrawer();
				break;
		}
		
		// clear existing diagram, if any
		var svg = document.getElementById(dst_id);
		while (svg.hasChildNodes()) 
		{
			svg.removeChild(svg.lastChild);
		}
		
		
		var cssStyle = document.createElementNS('http://www.w3.org/2000/svg','style');
		cssStyle.setAttribute('type','text/css');
		
		var style=document.createTextNode("text{font-size:6px;}");
		cssStyle.appendChild(style);
		
		svg.appendChild(cssStyle);
		
		
		var g = document.createElementNS('http://www.w3.org/2000/svg','g');
		g.setAttribute('id','viewport');
		svg.appendChild(g);
		
		
		td.Init(t, {svg_id: 'viewport', width:width, height:height, fontHeight:10, root_length:0.1} );
		
		td.CalcCoordinates();
		td.Draw();
		
		// label leaves...
		
		var n = new NodeIterator(t.root);
		var q = n.Begin();
		while (q != null)
		{
			if (q.IsLeaf())
			{
				switch (drawing_type)
				{
					case 'circle':
					case 'circlephylogram':
						var align = 'left';
						var angle = q.angle * 180.0/Math.PI;
						if ((q.angle > Math.PI/2.0) && (q.angle < 1.5 * Math.PI))
						{
							align = 'right';
							angle += 180.0;
						}
						document.getElementById('message').innerHTML='Drawing ' + q.label;
						drawRotatedText('viewport', q.xy, q.label, angle, align)
						break;
				
					case 'cladogram':
					case 'rectanglecladogram':
					case 'phylogram':
					default:				
						drawText('viewport', q.xy, q.label);
						document.getElementById('message').innerHTML='Drawing ' + q.label;
						break;
				}
			}
			q = n.Next();
		}
		
		document.getElementById('message').innerHTML='Done';
		// pan
		$("#" + dst_id).svgPan('viewport');
	}
}