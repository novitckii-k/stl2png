import struct
import time
import zlib

from basicgeo import P2, Partition1, along
from trianglebarmesh import TriangleBarMesh


class TriZSlice:
    def __init__(self, optionverbose):
        self.optionverbose = optionverbose
        self.tbms = []
        # multiple stlfiles will be unioned together
        # eg like support material structures,
        # though these may have a better definition than triangle files
        # such as branching lines and variable offset radii
        
    def load_stl_file(self, stlfile, transmap):
        # every edge has nodefrom.p.z<=nodeto.p.z
        tbm = TriangleBarMesh(stlfile, transmap, nodesortkey=lambda x: (x[0][2], x[0][1], x[0][0], x[1]))
        if self.optionverbose:
            nnodes, nedges, ntriangles, nsinglesidededges = tbm.get_facts()
            print("%s:  nodes=%d edges=%d triangles=%d singlesidededges=%d" % (stlfile, nnodes, nedges, ntriangles,
                                                                               nsinglesidededges))
            if nsinglesidededges != 0:
                print("*** Warning, not a closed surface")
        self.tbms.append(tbm)

    def set_extents(self, optionextra):
        self.xlo, self.xhi = min(tbm.xlo for tbm in self.tbms), max(tbm.xhi for tbm in self.tbms)
        self.ylo, self.yhi = min(tbm.ylo for tbm in self.tbms), max(tbm.yhi for tbm in self.tbms)
        self.zlo, self.zhi = min(tbm.zlo for tbm in self.tbms), max(tbm.zhi for tbm in self.tbms)
        if self.optionverbose:
            print("Dimensions: xlo=%.3f xhi=%.3f  ylo=%.3f yhi=%.3f  zlo=%.3f zhi=%.3f" % (self.xlo, self.xhi, self.ylo,
                                                                                           self.yhi, self.zlo,
                                                                                           self.zhi))
        if optionextra[-1] == "%":
            expc = float(optionextra[:-1])
            xex = (self.xhi-self.xlo)*expc/100
            yex = (self.yhi-self.ylo)*expc/100
        else:
            xex = float(optionextra)
            yex = float(optionextra)
        self.xlo -= xex;  self.xhi += xex
        self.ylo -= yex;  self.yhi += yex

    def build_pixel_grid_structures(self, optionwidthpixels, optionheightpixels):
        # grids partitions defining the interval width of each pixel
        self.xpixels = Partition1(self.xlo, self.xhi, optionwidthpixels)
        if optionheightpixels == 0:
            heightpixels = int(optionwidthpixels/(self.xhi - self.xlo)*(self.yhi - self.ylo) + 1)
            newyhi = heightpixels*(self.xhi - self.xlo)/optionwidthpixels + self.ylo
            if self.optionverbose:
                print("Revising yhi from %.3f to %.3f to for a whole number of pixels" % (self.yhi, newyhi))
            self.yhi = newyhi
        else:
            heightpixels = optionheightpixels
        self.ypixels = Partition1(self.ylo, self.yhi, heightpixels)

        xpixwid = self.xpixels.vs[1] - self.xpixels.vs[0]
        ypixwid = self.ypixels.vs[1] - self.ypixels.vs[0]
        
        if self.optionverbose:
            print("numxpixels=%d  each width=%.3f  x0=%0.3f" % (self.xpixels.nparts, xpixwid, self.xpixels.vs[0]))
            print("numypixels=%d  each height=%.3f  y0=%0.3f" % (self.ypixels.nparts, ypixwid, self.ypixels.vs[0]))
        
        # partitions with interval boundaries down middle of each pixel with extra line each side for convenience
        self.xpixmidsE = Partition1(self.xlo - xpixwid*0.5, self.xhi + xpixwid*0.5, self.xpixels.nparts + 1)
        self.ypixmidsE = Partition1(self.ylo - ypixwid*0.5, self.yhi + ypixwid*0.5, self.ypixels.nparts + 1)

    def calc_pixel_y_cuts(self, z, tbm):
        tbarpairs = []
        barcuts = {}
        for bar in tbm.bars:   # bucketing could speed this up
            assert bar.nodeback.p.z <= bar.nodefore.p.z
            if bar.nodeback.p.z <= z < bar.nodefore.p.z:
                bar1 = bar.barforeright
                node2 = bar1.get_node_fore(bar1.nodeback == bar.nodefore)
                bar_c = bar1 if node2.p.z <= z else bar1.get_fore_right_bl(bar1.nodeback == bar.nodefore)
                tbarpairs.append((bar.i, bar_c.i))
                
                lam = (z - bar.nodeback.p.z)/(bar.nodefore.p.z - bar.nodeback.p.z)
                cx = along(lam, bar.nodeback.p.x, bar.nodefore.p.x)
                cy = along(lam, bar.nodeback.p.y, bar.nodefore.p.y)
                barcuts[bar.i] = P2(cx, cy)

        ycuts = [[] for _ in range(self.ypixels.nparts)]  # plural for set of all raster rows
        for i, i1 in tbarpairs:
            p0, p1 = barcuts[i], barcuts[i1]
            iyl, iyh = self.ypixmidsE.get_part_range(min(p0.v, p1.v), max(p0.v, p1.v))
            for iy in range(iyl, iyh):
                yc = self.ypixmidsE.vs[iy+1]
                assert (p0.v < yc) != (p1.v < yc)
                lam = (yc - p0.v)/(p1.v - p0.v)
                xc = along(lam, p0.u, p1.u)
                ycuts[iy].append(xc)

        return ycuts
        
    def consolidate_y_cut_singular(self, ycutlist):
        l_ycuts = []
        for i, ycuts in enumerate(ycutlist):
            ycuts.sort()
            for j, yc in enumerate(ycuts):
                l_ycuts.append((yc, i, (j % 2 == 1)))
        l_ycuts.sort()
        l_i = set()
        ysegs = []
        yclo = -1.0
        for yc, i, bout in l_ycuts:
            if bout:
                l_i.remove(i)
                if len(l_i) == 0:
                    ychi = yc
                    ysegs.append((yclo, ychi))
            else:
                if len(l_i) == 0:
                    yclo = yc
                assert i not in l_i
                l_i.add(i)
        assert len(l_i) == 0
        return ysegs

    def calc_y_segrasters(self, z):
        ysegrasters = []
        ycuts_list = [self.calc_pixel_y_cuts(z, tbm) for tbm in self.tbms]
        for iy in range(self.ypixels.nparts):  # work through each raster line across the list of stlfiles
            ycutlist = [ycuts[iy] for ycuts in ycuts_list]
            ysegs = self.consolidate_y_cut_singular(ycutlist)
            ysegrasters.append(ysegs)
        return ysegrasters
        
    def calc_naked_compressed_bitmap(self, ysegrasters):
        compressor = zlib.compressobj()
        lcompressed = []

        def addcompressl(s):
            c = compressor.compress(s)
            if c:
                lcompressed.append(c)
        black_line = b"\0"*(self.xpixels.nparts+1)
        white_line = b"\xFF"*self.xpixels.nparts
        assert len(ysegrasters) == self.ypixels.nparts
        for ysegs in ysegrasters:
            previxhl = -1  # to get an extra \0 at the start
            for yseg in ysegs:
                ixl, ixh = self.xpixmidsE.get_part_range(yseg[0], yseg[1])
                addcompressl(black_line[:ixl-previxhl])
                addcompressl(white_line[:ixh-ixl])
                previxhl = ixh
            addcompressl(black_line[:self.xpixels.nparts-previxhl])
        lcompressed.append(compressor.flush())
        return lcompressed
            
    # this is a very low volume implementation of the PNG standard
    def write_png(self, fout, lcompressed):
        widthpixels, heightpixels = self.xpixels.nparts, self.ypixels.nparts
        fout.write(b"\x89" + "PNG\r\n\x1A\n".encode('ascii'))
        colortype, bitdepth, compression, filtertype, interlaced = 0, 8, 0, 0, 0
        block = struct.pack("!I4sIIBBBBB", 13, "IHDR".encode('ascii'), widthpixels, heightpixels, bitdepth, colortype,
                            compression, filtertype, interlaced)
        fout.write(block)
        fout.write(struct.pack("!I", zlib.crc32(block[4:]) & 0xFFFFFFFF))
        # length of compressed data at start of compressed section (which is why it can't be streamed out)
        fout.write(struct.pack("!I", sum(map(len, lcompressed))))
        idat = "IDAT".encode('ascii')
        crc = zlib.crc32(idat)
        fout.write(idat)
        for c in lcompressed:
            crc = zlib.crc32(c, crc)
            fout.write(c)
        fout.write(struct.pack("!I", crc & 0xFFFFFFFF))
        block = "IEND".encode('ascii')
        bcrc = zlib.crc32(block)
        fout.write(struct.pack("!I4sI", 0, block, bcrc & 0xFFFFFFFF))
        fout.close()

    def slice_to_png(self, z, pngname):
        stime = time.time()
        ysegrasters = self.calc_y_segrasters(z)
        lcompressed = self.calc_naked_compressed_bitmap(ysegrasters)
        self.write_png(open(pngname, "wb"), lcompressed)
        if self.optionverbose:
            print("Sliced at z=%f to file %s  compressbytes=%d %dms" % (z, pngname, sum(map(len, lcompressed)),
                                                                        (time.time()-stime)*1000))
        
        # conts = [ ]
        # for iy, ysegs in enumerate(ysegrasters):
        #    conts.extend([[(yseg[0], tzs.ypixmidsE.vs[iy+1]), (yseg[1], tzs.ypixmidsE.vs[iy+1])]  for yseg in ysegs])
        # sendactivity(contours=conts)
