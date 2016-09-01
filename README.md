# CLAS12-SVT-Geometry

<ul>
<li> SVTFactory
<ul>
<li> SVTConstants
<ul>
<li> Connects to CCDB and loads core parameters. </li>
<li> Option to load alignment shift data from file. </li>
<li> Provides conversions between indexing conventions. </li>
</ul>
</li>
<li> SVTStripFactory
<ul>
<li> Geometry factory for sensor strips. </li>
</ul>
</li>
<li> SVTVolumeFactory
<ul>
<li> Geometry factory for detector volumes. </li>
</ul>
</li>
<li> SVTAlignmentFactory
<ul>
<li> Geometry factory for fiducial points and file I/O for alignment data. </li>
</ul>
</li>
</ul>
</li>

<li> Alignment
<ul>
<li> AlignmentFactory
<ul>
<li> Universal class for processing and applying alignment shifts to points and volumes. </li>
</ul>
</li>
</ul>
</li>

<li>
Misc
<ul>
<li> Util
<ul>
<li> Universal utility class for vector and volume manipulation, rotation conversions, and file I/O. </li>
</ul>
</li>
<li> Matrix
<ul>
<li> Univerisal class for basic matrix algebra. </li>
<li> Supports conversions for 3D rotations: Tait-Bryan, axis-angle. </li>
</ul>
</li>
</ul>
</li>
</ul>
