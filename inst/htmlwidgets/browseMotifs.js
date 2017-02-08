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
        g = svg.append("g").attr("transform", "translate(40,0)");
        //zoom
        var zoom = d3.zoom()
                     .scaleExtent([1, 40])
                     .translateExtent([[-100, -100], 
                                      [width + 90, height + 100]])
                     .on("zoom", zoomed);
        function zoomed(){
          g.attr("transform", d3.event.transform);
        }
        svg.call(zoom);
        console.log(x);
        //leave top and bottom a full motif height
        //leave right+left 40+20+motif width
        var tree = d3[x.layout]()
                    .size([height-x.maxH, width - 60 - x.maxW]);
        //console.log(JSON.stringify(x));
        var root = d3.hierarchy(x.elements);
        //console.log(root);
        tree(root);
        
        var link = g.selectAll(".link")
                    .data(root.descendants().slice(1))
                    .enter().append("path")
                    .attr("class", "link")
                    .attr("d", function(d) {
                      return "M" + d.y + "," + d.x + 
                             "C" + (d.parent.y + 
                               Math.abs(d.parent.y-d.y)/2) + 
                             "," + d.x + " " + 
                             (d.parent.y + 
                               Math.abs(d.parent.y-d.y)/2) + 
                             "," + d.parent.x + 
                             " " + d.parent.y + 
                             "," + d.parent.x;
                    });

        var node = g.selectAll(".node")
              .data(root.descendants())
              .enter().append("g")
              .attr("class", function(d) {
                return "node" + (d.children ? 
                                  " node--internal" : 
                                  " node--leaf"); 
              })
              .attr("transform", function(d) {
                return "translate(" + d.y + "," + d.x + ")"; })
              .call(d3.drag()
                    .on("start", dragstarted)
                    .on("drag", dragged)
                    .on("end", dragended));
          function dragstarted(d) {
            d3.select(this).raise().classed("active", true);
          }
          function dragged(d) {
            d.x += d3.event.dy;
            d.y += d3.event.dx;
            var node = d3.select(this);
            node.attr("transform", function(d){
                return "translate(" + d.y + "," + d.x + ")";
            });
            var links = g.selectAll(".link");
            links.attr("d", function(d) {
                      return "M" + d.y + "," + d.x + 
                             "C" + (d.parent.y + 
                               Math.abs(d.parent.y-d.y)/2) + 
                             "," + d.x + " " + 
                             (d.parent.y + 
                               Math.abs(d.parent.y-d.y)/2) + 
                             "," + d.parent.x + 
                             " " + d.parent.y + 
                             "," + d.parent.x;
                    });
          }
          function dragended(d) {
            d3.select(this).classed("active", false);
          }

          console.log(node);
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
