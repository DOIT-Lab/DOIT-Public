NVAst=cell2struct(NVA(2:2:end), NVA(1:2:end), 2);

typ=split(typNm, '_');

switch typ{3}
    case 'I'
        datNm='\acrfull{I}';
    case 'P'
        datNm='\acrfull{phi}';
    case 'GI'
        datNm='gated \acrfull{I}';
    case 'T'
        datNm='\acrfull{mTOF}';
    case 'DGI'
        datNm='difference in gated \acrfull{I}';
    case 'V'
        datNm='\acrfull{var}';
    otherwise
end

fid=fopen(['Figs/' h.Name '.tex'], 'w+');

fprintf(fid, "\\begin{figure*}\n");
fprintf(fid, "\t\\begin{center}\n");
fprintf(fid, "\t\t\\includegraphics{\\gitPath makeFigs/Figs/%s_%d.pdf}\n", ...
    h.Name, h.Number);
fprintf(fid, "\t\\end{center}\n");
fprintf(fid, "\t\\caption{");
fprintf(fid, "$x$-$z$ plane of the \\acrfull{sen} to a ");
fprintf(fid, "$\\SI{%.1f}{\\milli\\meter}" + ...
    "\\times\\SI{%.1f}{\\milli\\meter}" + ...
    "\\times\\SI{%.1f}{\\milli\\meter}$ ", NVAst.pert);
fprintf(fid, "perturbation scanned \\SI{%.1f}{\\milli\\meter} ", NVAst.dr);
fprintf(fid, "measured by \\acrfull{%s} \\acrfull{%s} %s. ",...
    typ{1}, typ{2}, datNm);
fprintf(fid, "(a)-(%s) ", char(96+varNum));
fprintf(fid, "Different values of \\acrfull{%s}. ", varNm);
fprintf(fid, "Generated using \\acrfull{%s}.", NVAst.simTyp);
fprintf(fid, "\\\\ \n");
if ~strcmp(varNm, 'rho') && ...
    ~strcmp(varNm, 'mrho') && ...
    ~strcmp(varNm, 'drho')
    switch typ{2}
        case 'SD'
            fprintf(fid, "\t\\Acrfull{rho}: \\SI{%.1f}{\\milli\\meter}\\\\ \n",...
                vecnorm(rs-rd, 2, 2));
        case 'SS'
            fprintf(fid, "\t\\Acrfullpl{rho}: [%.1f, %.1f]" + ...
                "~\\si{\\milli\\meter}\\\\ \n",...
                vecnorm(rs-rd, 2, 2));
        case 'DS'
            fprintf(fid, "\t\\Acrfullpl{rho}: [%.1f, %.1f, %.1f, %.1f]" + ...
                "~\\si{\\milli\\meter}\\\\ \n",...
                vecnorm([rs; rs]-[rd([1, 1], :); rd([2, 2], :)], 2, 2));
        otherwise
    end
end
if strcmp(varNm, 'mrho')
    fprintf(fid, "\t\\Acrfull{drho}: \\SI{%.1f}{\\milli\\meter}\\\\ \n",...
        drho);
end
if strcmp(varNm, 'drho')
    fprintf(fid, "\t\\Acrfull{mrho}: \\SI{%.1f}{\\milli\\meter}\\\\ \n",...
        mrho);
end
fprintf(fid, "\t\\Acrfull{n} inside: \\num{%.3f};\\quad", optProp.nin);
fprintf(fid, "\t\\Acrfull{n} outside: \\num{%.3f}\\\\ \n", optProp.nout);
fprintf(fid, "\t\\Acrfull{musp}: \\SI{%.2f}{\\per\\milli\\meter};\\quad", ...
    optProp.musp);
if strcmp(NVAst.simTyp, 'MC')
    fprintf(fid, "\t\\Acrfull{g}: \\num{%.1f}\\\\ \n", optProp.g);
end
fprintf(fid, "\t\\Acrfull{mua}: \\SI{%.3f}{\\per\\milli\\meter}\\\\ \n", ...
    optProp.mua);
if strcmp(typ{1}, 'FD') && ~strcmp(varNm, 'fmod')
    fprintf(fid, "\t\\Acrfull{fmod}: \\SI{%.0f}{\\mega\\hertz}\\\\ \n", ...
        NVAst.fmod/1e6);
end
if strcmp(typ{3}, 'DGI') && ~strcmp(varNm, 'mtg') && ~strcmp(varNm, 'dtg')
    fprintf(fid, "\tEarly \\acrfull{t} gate: [%.0f, %.0f]" + ...
        "~\\si{\\pico\\second}\\\\ \n",...
        NVAst.tgE);
end
if (strcmp(typ{3}, 'GI') || strcmp(typ{3}, 'DGI')) && ~strcmp(varNm, 'mt') ...
    && ~strcmp(varNm, 'dt') && ~strcmp(varNm, 'mtg') && ~strcmp(varNm, 'dtg')
    fprintf(fid, "\t\\Acrfull{t} gate: [%.0f, %.0f]" + ...
        "~\\si{\\pico\\second}\\\\ \n",...
        NVAst.tg);
end
if strcmp(varNm, 'mt')
    fprintf(fid, "\t\\Acrfull{dt} gate: %.0f" + ...
        "~\\si{\\pico\\second}\\\\ \n",...
        dt);
end
if strcmp(varNm, 'dt')
    fprintf(fid, "\t\\Acrfull{mt} gate: %.0f" + ...
        "~\\si{\\pico\\second}\\\\ \n",...
        mt);
end
if strcmp(varNm, 'mtg')
    fprintf(fid, "\t\\Acrfull{dt} gates: %.0f" + ...
        "~\\si{\\pico\\second}\\quad",...
        dt);
    fprintf(fid, "\t\\Acrfull{dtg}: %.0f" + ...
        "~\\si{\\pico\\second}\\\\ \n",...
        dtg);
end
if strcmp(varNm, 'dtg')
    fprintf(fid, "\t\\Acrfull{dt} gates: %.0f" + ...
        "~\\si{\\pico\\second}\\quad",...
        dt);
    fprintf(fid, "\t\\Acrfull{mtg}: %.0f" + ...
        "~\\si{\\pico\\second}\\\\ \n",...
        mtg);
end
if strcmp(NVAst.simTyp, 'MC')
    fprintf(fid, "\tDetector Numerical Aperature (NA): " + ...
        "\\num{%.1f};\\quad", NVAst.detNA);
    fprintf(fid, "\tNumber of photons: " + ...
        "\\num{%d}\\\\ \n", NVAst.np);
end
fprintf(fid, "\t}\\label{fig:%s}\n", h.Name);
fprintf(fid, "\\end{figure*}\n");

fopen(['Figs/' h.Name '.tex']);