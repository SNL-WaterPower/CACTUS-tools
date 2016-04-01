function T=CreateTurbineAndWriteGeom(GeomFN,varargin)

        T = CreateTurbine(varargin{:});
        WriteTurbineGeom(GeomFN, T)

end