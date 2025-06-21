class Chart{
  constructor(ctx,config){
    this.ctx=ctx;
    this.config=config;
    this.options=config.options||{};
    this.updateScales();
    this.draw();
  }
  updateScales(){
    const datasets=this.config.data.datasets;
    const xs=[],ys=[];
    datasets.forEach(ds=>{
      (ds.data||[]).forEach(p=>{xs.push(p.x); ys.push(p.y);});
    });
    this.minX=Math.min(...xs); this.maxX=Math.max(...xs);
    this.minY=Math.min(...ys); this.maxY=Math.max(...ys);
    const sx=x=>(x-this.minX)/(this.maxX-this.minX||1)*this.ctx.canvas.width;
    const sy=y=>(this.maxY-y)/(this.maxY-this.minY||1)*this.ctx.canvas.height;
    this.scales={x:{getPixelForValue:sx},y:{getPixelForValue:sy}};
  }
  destroy(){
    this.ctx.clearRect(0,0,this.ctx.canvas.width,this.ctx.canvas.height);
  }
  draw(){
    const ctx=this.ctx;
    ctx.clearRect(0,0,ctx.canvas.width,ctx.canvas.height);
    const sx=this.scales.x.getPixelForValue;
    const sy=this.scales.y.getPixelForValue;
    (this.config.data.datasets||[]).forEach(ds=>{
      ctx.strokeStyle=ds.borderColor||'black';
      ctx.fillStyle=ds.backgroundColor||'transparent';
      const data=ds.data||[];
      if(ds.type==='scatter'){
        data.forEach(p=>{ctx.beginPath();ctx.arc(sx(p.x),sy(p.y),3,0,2*Math.PI);ctx.fill();});
      }else{
        ctx.beginPath();
        data.forEach((p,i)=>{i?ctx.lineTo(sx(p.x),sy(p.y)):ctx.moveTo(sx(p.x),sy(p.y));});
        if(ds.showLine!==false) ctx.stroke();
        if(ds.fill){
          ctx.lineTo(sx(data[data.length-1].x),sy(0));
          ctx.lineTo(sx(data[0].x),sy(0));
          ctx.fill();
        }
      }
    });
    (Chart.plugins||[]).forEach(p=>p.afterDatasetsDraw&&p.afterDatasetsDraw(this,{},{}));
  }
}
Chart.plugins=[];
Chart.register=function(p){Chart.plugins.push(p);};
