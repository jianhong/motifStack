/*
 plotWebLogo
Copyright 2016, Jianhong Ou
E-mail: jianhong.ou@umassmed.edu
Licensed under the MIT license
*/
;(function(){
  d3.selection.prototype.plotWebLogo = function(options){
    var defaults = {
      x:0,
      y:0,
      height:100,
      width:300,
      margin:{bottom:5, left:5, top:5, right:5},
      letters:["A","C","G","T"],
      //color
      col:{"A":"#008000", "C":"#00c", "G":"#ffb300", "T":"#c00"},
      //background
      p:{"A":0.25, "C":0.25, "G":0.25, "T":0.25},
      xlab:"position",
      xaxis: true,
      ylab:"bits",
      yaxis: true,
      title:"jWebLogo",
      names:[],
      fontFamily: "Helvetica, Arial, sans-serif",
      fontSize: "12pt",
      defaultColor: "black",
      background:"transparent"
    };
    var self = this;
    if(typeof(options)=="undefined") return null;
    if(!options.hasOwnProperty("motif")) return null;
    if(!d3) return null;
    for(var key in defaults){
      if(!options.hasOwnProperty(key)){
        options[key] = defaults[key];
      }
    }
    for(var i=0; i<options.letters.length; i++){
          if(!options.col.hasOwnProperty(options.letters[i])) return null;
          if(!options.p.hasOwnProperty(options.letters[i])) return null;
    }
    //if(options.height<50) options.height=50;
    if(Object.prototype.toString.call( options.names ) !== '[object Array]'){
      options.names = [options.names];
      len = options.motif.length;
      while(options.names.length * 2 < len) options.names = options.names.concat(options.names);
      if(options.names.length < len) options.names = options.names.concat(options.names.slice(0, len - options.names.length));
    }
    
      function log2(x){
            if(x===0) x=1/1000000;
            return Math.log(x)/Math.log(2);
      }
    var M=options.motif;
//    console.log(M);
    var totalMotif=M.length;
    var totalStack=0;
    for(var i=0; i<totalMotif; i++){
      if(totalStack<M[i].length) totalStack=M[i].length;
    }
    var entropy = parseInt(10*log2(options.letters.length));
    var motifNameSpace=24;
    var defaultHeight=100;
    var titleHeight=options.title.length>0?24:0;
    var bottomHeight=options.xlab.length>0?24:0;
    var xtickHeight=options.xaxis?10:0;
    var totalHeight=(defaultHeight+motifNameSpace+xtickHeight)*totalMotif+
                     titleHeight+bottomHeight+options.margin.top+options.margin.bottom;
    //vertical align center
    options.y = options.y - totalHeight/2;
    //horizontal align left, gap=5
    options.x = options.x + 5;
    var ystop=options.margin.top+defaultHeight+motifNameSpace+xtickHeight;
    var yintv=defaultHeight/entropy;
    var zero_position={"x":options.margin.left,"y":ystop};
    var curr_position={"x":options.margin.left,"y":ystop};
    var xintv=parseInt(defaultHeight*0.6);
    var halfUnit=xintv/2;
    var paddingLeft=options.yaxis?xintv:0;
    var width=options.margin.left+paddingLeft+totalStack*xintv+options.margin.right;
    var ratios=[options.width/width, options.height/totalHeight];
    //color
    var black="#000";
    var white="#fff";
    self.each(function(){
      var $this = d3.select(this);
      var paper =  $this.append("g");
      //console.log(paper_pos);
      function drawBorder(){
        var rect= paper.append("rect")
                       .attr("x", 0)
                       .attr("y", 0)
                       .attr("width", width)
                       .attr("height",totalHeight)
                       .style("fill", options.background)
                       .style("color", options.defaultColor);
      }
      
      function line(x1,y1,x2,y2,color=options.defaultColor){
        return paper.append("line")
                    .attr("x1", x1)
                    .attr("y1", y1)
                    .attr("x2", x2)
                    .attr("y2", y2)
                    .style("fill", color)
                    .style("stroke", color);
      }
      
      function text(x, y, label, color=options.defaultColor, fontSize=options.fontSize, fontWeight="normal"){
        return paper.append("text")
                    .attr("x", x)
                    .attr("y", y)
                    .text(label)
                    .style("font-family", options.fontFamily)
                    .style("font-size", fontSize)
                    .style("font-weight", fontWeight)
                    .style("fill", color)
                    .style("text-anchor", "middle")
                    .style("dominant-baseline", "central");
      }
      
      function drawYaxis(){
        //draw axis
        var yaxis = line(zero_position.x, zero_position.y-defaultHeight,
                zero_position.x, zero_position.y);
        
        //draw ticks
        for(var i=0; i<entropy+1; i++){
          xlen=5;      
          if(i%10===0){
            xlen=10;
            var yticklab=text(zero_position.x-xlen-14, 
                              curr_position.y, 
                              i/10+".0",
                              "black",
                              "10pt");
          }
          var ytick=line(zero_position.x-xlen, curr_position.y, 
                         zero_position.x, curr_position.y);
          curr_position.y-=yintv;
        }
        //draw label
        var ylab=text(curr_position.x-48, curr_position.y+yintv*entropy/2, options.ylab);
        ylab.attr("transform", "rotate(-90 "+(curr_position.x-48)+" "+(curr_position.y+yintv*entropy/2)+")");
      }
      //object site
      function isInt(n){
        return typeof n === 'number' && parseFloat(n) == parseInt(n, 10) && !isNaN(n);
      }
      function SITE(x){
            function add(a, b){ return a + b;}
            var total=x.reduce(add, 0);
            var s = [];
            if(total==1){//is pfm
                  for(var i=0; i<options.letters.length; i++){
                        s[options.letters[i]] = x[i];
                  }
            }else{
                  for(var i=0; i<options.letters.length; i++){
                        s[options.letters[i]] = x[i]/total;
                  }
            }
            return(s);
      }
      function getIC(s){
        //calculate IC
        //formular: Ii=sum(f*log(f/p,base=2)['A,C,G,T'])
        var IC=0;
        for(var k in s){
          IC+=s[k]*(log2(s[k])-log2(options.p[k]));
        }
        return IC;
      }
      function stack(s){
        var ic=getIC(s);
        ic1=defaultHeight*ic/log2(options.letters.length);
        //calculate height
        var h=[];
        for(var k in s){
          h.push([k,s[k]*ic1,5*s[k]*ic/6]);
        }
        h.sort(function(a,b){return a[1] - b[1]});
        return h;
      }
      function drawStack(num,h){
            var x = curr_position.x;
        if(options.xaxis){
          //draw num, underline
          if(num!==0){
            var stacklab=text(x+halfUnit, zero_position.y+xtickHeight, num,
                            "black", "10pt");
          }
          //draw beam
          var stackline = line(x+2, zero_position.y, 
                               x+xintv-2, zero_position.y);
        }
        //draw stack
        if(num!==0){
          var pos_y=zero_position.y;
  //        console.log(h);
          for(var i=0; i<h.length; i++){
            if(h[i][1]!==0){
              var stack=text(halfUnit, 0, h[i][0],
                             options.col[h[i][0]], xintv+"pt", "bold");
              stack.attr("transform", 
                         "translate("+x+","+
                                      (pos_y-parseInt(h[i][1]/2))+
                         ") scale("+1+","+h[i][2]+")");
              pos_y-=h[i][1];
            }
          }
        }
        curr_position.x+=xintv;
      }
      function drawXlab(){
        //draw x label
        var xlab=text((width-zero_position.x)/2+zero_position.x,totalHeight-bottomHeight/2,options.xlab);
      }
      function drawTitle(){
        var title=text(width/2,titleHeight/2,options.title, "black", "14pt");
      }
      function drawMotifName(mn){
        var motifname=text((width-zero_position.x)/2+zero_position.x,zero_position.y-defaultHeight-motifNameSpace/2,mn,
                            "black", "10pt");
      }
      drawBorder();
      drawTitle();
      for(var i=0; i<totalMotif; i++){
        var motif=M[i];
        //renew zero_position;
        if(options.yaxis){
              zero_position.x=options.margin.left+paddingLeft;
                  zero_position.y=options.margin.top+titleHeight+(i+1)*(defaultHeight+motifNameSpace)+i*xtickHeight;
                  curr_position=zero_position;
                  drawYaxis();
        }
        
        //renew zero_position;
        zero_position.x=options.margin.left+paddingLeft;
        zero_position.y=options.margin.top+titleHeight+(i+1)*(defaultHeight+motifNameSpace)+i*xtickHeight;
        curr_position=zero_position;
        drawXlab();
        drawMotifName(options.names[i]);
        //draw stack
        //console.log(curr_position.x);
        for(var j=0; j<motif.length; j++){
          var s=new SITE(motif[j]);
          s=stack(s);
          drawStack(j+1,s);
          //console.log(curr_position.x);
        }
        if(motif.length<totalStack){
          for(var j=motif.length; j<totalStack; j++){
            drawStack(0,null);
          }
        }
      }
      paper.attr("transform", "scale("+ratios[0]+","+ratios[1]+") translate("+options.x+","+options.y+")");
    });
    return self;
  };
})();