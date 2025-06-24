const { useState, useEffect, useRef } = React;

function LineChart({xs, ys}) {
  const ref = useRef(null);
  useEffect(() => {
    const svg = d3.select(ref.current);
    svg.selectAll('*').remove();
    const width = 300;
    const height = 150;
    svg.attr('viewBox', `0 0 ${width} ${height}`);
    const x = d3.scaleLinear().domain(d3.extent(xs)).range([40, width - 10]);
    const y = d3.scaleLinear().domain(d3.extent(ys)).range([height - 20, 10]);
    const line = d3.line().x((d,i) => x(xs[i])).y(d => y(d));
    svg.append('path').attr('fill','none').attr('stroke','steelblue').attr('d', line(ys));
    svg.append('g').attr('transform',`translate(0,${height-20})`).call(d3.axisBottom(x));
    svg.append('g').attr('transform','translate(40,0)').call(d3.axisLeft(y));
  }, [xs, ys]);
  return React.createElement('svg', {ref, className:'w-full h-40 bg-white rounded'});
}

function BeamCalculator() {
  const [spans, setSpans] = useState([{length:5}]);
  const [results, setResults] = useState(null);

  function addSpan(){ setSpans([...spans,{length:5}]); }
  function updateSpan(idx,val){ const arr=[...spans]; arr[idx].length=parseFloat(val); setSpans(arr); }
  function removeSpan(idx){ setSpans(spans.filter((_,i)=>i!==idx)); }

  useEffect(()=>{
    const state={spans:spans.map(s=>({length:parseFloat(s.length)||0})),pointLoads:[],lineLoads:[],maxNodeDist:0.2,E:210e9,I:1e-6};
    const data=window.computeResults(state); if(!data) return;
    const diags=window.computeDiagrams(state,data.nodes,data.reactions);
    setResults({data,diags});
  },[spans]);

  function renderSpan(s,idx){
    return React.createElement('div',{key:idx,id:`span-${idx}`,className:'flex items-center space-x-2'},[
      React.createElement('label',{className:'w-32'},`Span ${idx+1} length`),
      React.createElement('input',{type:'number',value:s.length,min:0.1,step:0.1,onChange:e=>updateSpan(idx,e.target.value),className:'border px-2 w-24'}),
      React.createElement('button',{onClick:()=>removeSpan(idx),className:'text-red-500'},'Remove')
    ]);
  }

  return React.createElement('div',{className:'p-4',id:'container'},[
    React.createElement('div',{id:'controls',className:'space-y-4'},[
      React.createElement('h1',{className:'text-xl font-bold'},'Beam Calculator'),
      ...spans.map(renderSpan),
      React.createElement('button',{id:'addSpanBtn',onClick:addSpan,className:'px-2 py-1 bg-blue-500 text-white rounded'},'Add Span'),
      React.createElement('button',{id:'exportBtn',onClick:()=>navigator.clipboard.writeText(JSON.stringify({spans})),className:'px-2 py-1 bg-green-500 text-white rounded'},'Export')
    ]),
    results && React.createElement('div',{id:'charts',className:'flex-1 grid grid-cols-1 gap-4 mt-4'},[
      React.createElement('div',{className:'section'},React.createElement(LineChart,{xs:results.diags.xs,ys:results.diags.shearVals})),
      React.createElement('div',{className:'section'},React.createElement(LineChart,{xs:results.diags.xs,ys:results.diags.momentVals}))
    ])
  ]);
}

ReactDOM.createRoot(document.getElementById('root')).render(React.createElement(BeamCalculator));
