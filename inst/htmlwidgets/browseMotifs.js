HTMLWidgets.widget({

  name: 'browseMotifs',

  type: 'output',

  factory: function(el, width, height) {
    var svg = d3.select(el).append("svg");
    var menu = d3.select(el).append("input")
                 .attr("type", "button")
                 .attr("name", "export")
                 .attr("id", "export")
                 .attr("value", "exportSVG");
    return {
      renderValue: function(x) {
        function writeDownloadLink(){
            function fireEvent(obj,evt){
              var fireOnThis = obj;
              var evObj;
              if( document.createEvent ) {
                evObj = document.createEvent('MouseEvents');
                evObj.initEvent( evt, true, false );
                fireOnThis.dispatchEvent( evObj );
              } else if( document.createEventObject ) {
                evObj = document.createEventObject();
                fireOnThis.fireEvent( 'on' + evt, evObj );
              }
            }
            svgAsDataUri(svg.node(), 'motifStack.svg', function(uri){
                var a = document.createElement('a');
                a.href = uri;
                a.download = 'motifStack.svg';
                fireEvent(a, 'click');
            });
        }
        d3.select("#export")
          .on("click", writeDownloadLink);
        
        svg.attr("width", width)
           .attr("height", height);
        //leave left side 40px space for root name
        
        var g = svg.append("g");
        if(x.layout!="radialPhylog"){
          g.attr("transform", "translate(40,0)");
        }
        //zoom
        var zoom = d3.zoom()
                     .scaleExtent([0, 40])
                     .on("zoom", zoomed);
        function zoomed(){
          g.attr("transform", d3.event.transform);
        }
        svg.call(zoom);
        //console.log(x);
        //leave top and bottom a full motif height
        //leave right+left 40+20+motif width
        var tree;
        var rw = width/2 - x.maxW, rh = height/2 - x.maxW;
        if(x.layout==="radialPhylog"){
          if(rw>rh){
            var tmp = rw;
            rw = rw/1.6;
            rh = tmp;
          }else{
            rw = rh/1.6;
          }
          tree = d3.cluster()
                   .size([rw, rh])
                   .separation(function(a, b) { 
                     return (a.parent == b.parent ? 1 : 2) / a.depth; 
                   });
        }else{
          tree = d3[x.layout]()
                    .size([height-x.maxH, width - 60 - x.maxW]);
        }
        //console.log(JSON.stringify(x));
        var root = d3.hierarchy(x.elements);
        //console.log(root);
        tree(root);
        function cartesian2Polar(x, y, dx, dy){
          var xy = project(x, y);
          x = xy[0] - width/2 + dx;
          y = xy[1] - height/2 + dy;
          radius = Math.sqrt(x*x + y*y);
          angle = Math.atan2(y,x); //This takes y first
          y = radius;
          x = angle * rw / 2 / Math.PI + rw/4;
          return [x, y];
        }
        function project(x, y) {
          var angle = (x - (rw/4))/ (rw/2) * Math.PI, radius = y;
          return [radius * Math.cos(angle)+width / 2, radius * Math.sin(angle)+height / 2];
        }
        function diagonal(y0, y1){
          return y0+Math.abs(y1-y0)/2;
        }
        var linkFun = {"cluster":function(d){
                                  return "M" + d.y + "," + d.x + 
                                         "C" + diagonal(d.parent.y, d.y) + 
                                         "," + d.x + " " + 
                                         diagonal(d.parent.y, d.y) + 
                                         "," + d.parent.x + 
                                         " " + d.parent.y + 
                                         "," + d.parent.x;
                                }, 
                       "tree": function(d){
                                  return "M" + d.y + "," + d.x + 
                                         "C" + diagonal(d.parent.y, d.y) + 
                                         "," + d.x + " " + 
                                         diagonal(d.parent.y, d.y) + 
                                         "," + d.parent.x + 
                                         " " + d.parent.y + 
                                         "," + d.parent.x;
                                }, 
                       "radialPhylog":function(d){
                                      return "M" + project(d.x, d.y) + 
                                             "C" + project(d.x, (d.y + d.parent.y) / 2) + 
                                             " " + project(d.parent.x, (d.y + d.parent.y) / 2) + 
                                             " " + project(d.parent.x, d.parent.y);
                       }};
        var link = g.selectAll(".link")
                    .data(root.descendants().slice(1))
                    .enter().append("path")
                    .attr("class", "link")
                    .attr("d", linkFun[x.layout]);
        var nodeFun = {
          "cluster":function(d) {
                return "translate(" + d.y + "," + d.x + ")"; },
          "tree":function(d) {
                return "translate(" + d.y + "," + d.x + ")"; },
          "radialPhylog":function(d) {
                return "translate(" + project(d.x, d.y) + ") rotate(" + (d.x < width/2 ? d.x*360/rw - 90 : d.x*360/rw + 90) + ")"; }
        };
        var node = g.selectAll(".node")
              .data(root.descendants())
              .enter().append("g")
              .attr("class", function(d) {
                return "node" + (d.children ? 
                                  " node--internal" : 
                                  " node--leaf"); 
              })
              .attr("transform", nodeFun[x.layout])
              .call(d3.drag()
                    .on("start", dragstarted)
                    .on("drag", dragged)
                    .on("end", dragended));
          function dragstarted(d) {
            d3.select(this).raise().classed("active", true);
          }
          function dragged(d) {
            if(x.layout==="radialPhylog"){
              var dxy = cartesian2Polar(d.x, d.y, d3.event.dx, d3.event.dy);
              d.x = dxy[0];
              d.y = dxy[1];
            }else{
              d.x += d3.event.dy;
              d.y += d3.event.dx;
            }
            var node = d3.select(this);
            node.attr("transform", nodeFun[x.layout]);
            var links = g.selectAll(".link");
            links.attr("d", linkFun[x.layout]);
          }
          function dragended(d) {
            d3.select(this).classed("active", false);
          }

          //console.log(node);
          node.append("circle").attr("r", x.nodeRadius);
        
          g.selectAll(".node--leaf").append("text")
              .attr("dy", 3)
              .attr("x", -8)
              .style("text-anchor", "end")
              .text(function(d) { return d.data.name; });
        
        node.each(function(d){
            var canvas=d3.select(this);
            if(!d.children){
                canvas.plotWebLogo({motif: d.data.motif,
                            height: d.data.height,
                            width: d.data.width,
                            names: "",
                            title: "",
                            xlab: "",
                            letters: d.data.letters,
                            col:d.data.color,
                              //background
                              p:d.data.background,
                            xaxis:d.data.xaxis,
                            yaxis:d.data.yaxis,
                            x: 0,
                            y: 0});
            }
        });
      },

      resize: function(width, height) {
        svg.attr("width", width)
           .attr("height", height);
      },

      svg: svg
    };
  }
});
