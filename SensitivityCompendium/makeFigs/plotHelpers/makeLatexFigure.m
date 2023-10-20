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

if strcmp(plotPreset, 'vox')
    fid=fopen(['Figs/' h.Name '_vox.tex'], 'w+');
else
    fid=fopen(['Figs/' h.Name '.tex'], 'w+');
end

fprintf(fid, "\\begin{figure*}\n");
fprintf(fid, "\t\\begin{center}\n");
fprintf(fid, "\t\t\\includegraphics{\\gitPath makeFigs/Figs/%s_%d.pdf}\n", ...
    h.Name, h.Number);
fprintf(fid, "\t\\end{center}\n");
fprintf(fid, "\t\\caption{");
fprintf(fid, "Third angle projection of the \\acrfull{sen} to a ");
switch plotPreset
    case 'none'
        fprintf(fid, "$\\SI{%.1f}{\\milli\\meter}" + ...
            "\\times\\SI{%.1f}{\\milli\\meter}" + ...
            "\\times\\SI{%.1f}{\\milli\\meter}$ ", NVAst.pert);
    case 'vox'
        fprintf(fid, "$\\SI{%.1f}{\\milli\\meter}" + ...
            "\\times\\SI{%.1f}{\\milli\\meter}" + ...
            "\\times\\SI{%.1f}{\\milli\\meter}$ ", ...
            NVAst.dr, NVAst.dr, NVAst.dr);
    otherwise
end
fprintf(fid, "perturbation scanned \\SI{%.1f}{\\milli\\meter} ", NVAst.dr);
fprintf(fid, "measured by \\acrfull{%s} \\acrfull{%s} %s. ",...
    typ{1}, typ{2}, datNm);
fprintf(fid, "(a) $x$-$y$ plane sliced at $z=\\SI{%.1f}{\\milli\\meter}$. ",...
    zsl);
fprintf(fid, "(b) Iso-surface sliced at ");
switch plotPreset
    case 'none'
        fprintf(fid, "$\\as{S}=%.3f$", posSurf);
        if ~isnan(negSurf)
            fprintf(fid, " and $\\as{S}=%.3f$. ", negSurf);
        else
            fprintf(fid, ". ");
        end
    case 'vox'
        fprintf(fid, "$\\as{S}=%.3f\\times 10^{-5}$", posSurf*1e5);
        if ~isnan(negSurf)
            fprintf(fid, "and $\\as{S}=%.3f\\times 10^{-5}$. ", negSurf*1e5);
        else
            fprintf(fid, ". ");
        end
    otherwise
end
fprintf(fid, "(c) $x$-$z$ plane sliced at $y=\\SI{%.1f}{\\milli\\meter}$. ",...
    ysl);
fprintf(fid, "(d) $y$-$z$ plane sliced at $x=\\SI{%.1f}{\\milli\\meter}$. ",...
    xsl);
fprintf(fid, "Generated using \\acrfull{%s}.", NVAst.simTyp);
fprintf(fid, "\\\\ \n");
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
fprintf(fid, "\t\\Acrfull{n} inside: \\num{%.3f};\\quad", optProp.nin);
fprintf(fid, "\t\\Acrfull{n} outside: \\num{%.3f}\\\\ \n", optProp.nout);
fprintf(fid, "\t\\Acrfull{musp}: \\SI{%.2f}{\\per\\milli\\meter};\\quad", ...
    optProp.musp);
if strcmp(NVAst.simTyp, 'MC')
    fprintf(fid, "\t\\Acrfull{g}: \\num{%.1f}\\\\ \n", optProp.g);
end
fprintf(fid, "\t\\Acrfull{mua}: \\SI{%.3f}{\\per\\milli\\meter}\\\\ \n", ...
    optProp.mua);
if strcmp(typ{1}, 'FD')
    fprintf(fid, "\t\\Acrfull{fmod}: \\SI{%.0f}{\\mega\\hertz}\\\\ \n", ...
        NVAst.fmod/1e6);
end
if strcmp(typ{3}, 'DGI')
    fprintf(fid, "\tEarly \\acrfull{t} gate: [%.0f, %.0f]" + ...
        "~\\si{\\pico\\second};\\quad\n",...
        NVAst.tgE);
end
if strcmp(typ{3}, 'GI') || strcmp(typ{3}, 'DGI')
    fprintf(fid, "\t\\Acrfull{t} gate: [%.0f, %.0f]" + ...
            "~\\si{\\pico\\second}\\\\ \n",...
            NVAst.tg);
end
if strcmp(NVAst.simTyp, 'MC')
    fprintf(fid, "\tDetector Numerical Aperature (NA): " + ...
        "\\num{%.1f};\\quad", NVAst.detNA);
    fprintf(fid, "\tNumber of photons: " + ...
        "\\num{%d}\\\\ \n", NVAst.np);
end
if strcmp(plotPreset, 'vox')
    fprintf(fid, "\t}\\label{fig:%s_vox}\n", h.Name);
else
    fprintf(fid, "\t}\\label{fig:%s}\n", h.Name);
end
fprintf(fid, "\\end{figure*}\n");

fopen(['Figs/' h.Name '.tex']);