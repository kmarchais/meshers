import * as THREE from "three";

// Conservative field stepping, not a signed-distance function. Rendering only:
// pixel tolerance and the finite step budget are not geometric error guarantees.
export function createPreview(config, cut) {
  const material = new THREE.ShaderMaterial({
    side: THREE.BackSide,
    transparent: true,
    uniforms: {
      shape: { value: config.shape },
      repeats: { value: config.repeat },
      thickness: { value: config.thickness },
      grade: { value: config.grade },
      cut: { value: cut },
    },
    vertexShader: `varying vec3 positionWorld;
      void main(){positionWorld=(modelMatrix*vec4(position,1.)).xyz;
      gl_Position=projectionMatrix*viewMatrix*vec4(positionWorld,1.);}`,
    fragmentShader: `
      varying vec3 positionWorld;
      uniform int shape;
      uniform float repeats, thickness, grade, cut;
      float band(vec3 p){
        vec3 q=p*repeats*6.28318530718;
        vec3 s=sin(q),c=cos(q),s2=sin(2.*q),c2=cos(2.*q);
        float f;
        if(shape==0) f=dot(s,c.yzx);
        else if(shape==1) f=1.1*dot(s2,c.yzx*s.zxy)
          -.2*dot(c2,c2.yzx)-.4*(c2.x+c2.y+c2.z);
        else f=c.x+c.y+c.z;
        return abs(f)-.5*thickness*(1.+2.*grade*p.x);
      }
      vec4 trace(vec3 target){
        vec3 ro=cameraPosition,rd=normalize(target-ro);
        vec3 low=vec3(-.5),high=vec3(min(.5,cut),.5,.5);
        if(high.x<=low.x) return vec4(0.);
        // Protect exactly parallel rays without changing their meaningful sign.
        vec3 safe=vec3(abs(rd.x)<1e-8?1e-8:rd.x,abs(rd.y)<1e-8?1e-8:rd.y,abs(rd.z)<1e-8?1e-8:rd.z);
        vec3 a=(low-ro)/safe,b=(high-ro)/safe;
        vec3 nearT=min(a,b),farT=max(a,b);
        float entry=max(max(nearT.x,nearT.y),nearT.z);
        float end=min(min(farT.x,farT.y),farT.z);
        float t=max(entry,0.);
        if(end<t) return vec4(0.);
        vec3 p=ro+rd*t;
        bool cap=entry>=0. && band(p)<0.;
        bool inside=band(p)<0.;
        bool hit=cap;
        float L=(shape==0?22.:shape==1?140.:11.)*repeats+thickness*abs(grade);
        for(int i=0;i<768;i++){
          if(hit||t>end) break;
          p=ro+rd*t;
          float f=band(p);
          if(abs(f)<.0004 || ((f<0.)!=inside)){hit=true;break;}
          t+=max(abs(f)/L,.000025);
        }
        if(!hit||t>end) return vec4(0.);
        vec3 n;
        if(cap){
          vec3 d=min(abs(p-low),abs(p-high));
          if(d.x<=d.y&&d.x<=d.z)n=vec3(abs(p.x-low.x)<abs(p.x-high.x)?-1.:1.,0.,0.);
          else if(d.y<=d.z)n=vec3(0.,abs(p.y-low.y)<abs(p.y-high.y)?-1.:1.,0.);
          else n=vec3(0.,0.,abs(p.z-low.z)<abs(p.z-high.z)?-1.:1.);
        }else{
          float h=.0001/repeats;
          n=normalize(vec3(band(p+vec3(h,0,0))-band(p-vec3(h,0,0)),band(p+vec3(0,h,0))-band(p-vec3(0,h,0)),band(p+vec3(0,0,h))-band(p-vec3(0,0,h))));
        }
        if(dot(n,rd)>0.)n=-n;
        float light=.27+.64*max(dot(n,normalize(vec3(2.,4.,3.))),0.)+.2*max(dot(n,normalize(vec3(-3.,1.,-2.))),0.);
        vec3 color=vec3(.69,.82,.51)*light;
        return vec4(color,1.);
      }
      void main(){
        // Four subpixel rays smooth implicit edges that raster MSAA cannot see.
        vec3 dx=dFdx(positionWorld)*.25,dy=dFdy(positionWorld)*.25;
        vec4 sum=trace(positionWorld-dx-dy)+trace(positionWorld+dx-dy)
                +trace(positionWorld-dx+dy)+trace(positionWorld+dx+dy);
        if(sum.a==0.) discard;
        gl_FragColor=vec4(sum.rgb/sum.a,sum.a*.25);
        #include <tonemapping_fragment>
        #include <colorspace_fragment>
      }`,
  });
  return new THREE.Mesh(new THREE.BoxGeometry(1, 1, 1), material);
}
